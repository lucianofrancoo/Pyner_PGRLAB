from __future__ import annotations

import json
import re
import unicodedata
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple


VALID_LABELS = {"alta", "media", "baja"}
VALID_EVIDENCE = {"directa", "indirecta", "débil"}


@dataclass
class ClassificationContext:
    organisms: List[str]
    strategies: List[str]
    conditions: List[str]
    tissues: List[str]
    free_terms: List[str]


class ResultClassifier:
    """Clasifica resultados de PubMed/BioProject usando Ollama o heuristica."""

    version = "1.0.0"

    def __init__(self, ollama_client: Optional[Any] = None):
        self.ollama_client = ollama_client

    @property
    def llm_available(self) -> bool:
        if not self.ollama_client:
            return False
        try:
            return bool(self.ollama_client.is_available())
        except Exception:
            return False

    def classify_pubmed(self, record: Dict[str, Any], query_payload: Dict[str, Any], use_llm: bool) -> Dict[str, Any]:
        context = self._build_context(query_payload)
        heuristic = self._heuristic_pubmed(record, context)

        if use_llm and self.llm_available:
            llm = self._classify_with_ollama(source="pubmed", record=record, context=context)
            if llm:
                return {**record, "classification": llm}

        return {**record, "classification": heuristic}

    def classify_bioproject(self, record: Dict[str, Any], query_payload: Dict[str, Any], use_llm: bool) -> Dict[str, Any]:
        context = self._build_context(query_payload)
        heuristic = self._heuristic_bioproject(record, context)

        if use_llm and self.llm_available:
            llm = self._classify_with_ollama(source="bioproject", record=record, context=context)
            if llm:
                return {**record, "classification": llm}

        return {**record, "classification": heuristic}

    def _build_context(self, payload: Dict[str, Any]) -> ClassificationContext:
        extracted = payload.get("extracted", {}) or {}
        synonyms = payload.get("synonyms", {}) or {}

        organisms = self._unique(
            [
                extracted.get("organism"),
                *(extracted.get("organism_variants", []) or []),
                *(synonyms.get("organism", []) or []),
            ]
        )
        strategies = self._unique(
            [
                *(extracted.get("strategies", []) or []),
                *(synonyms.get("strategies", []) or []),
            ]
        )
        conditions = self._unique(
            [
                *(extracted.get("conditions", []) or []),
                *(synonyms.get("conditions", []) or []),
            ]
        )
        tissues = self._unique(
            [
                *(extracted.get("tissues", []) or []),
                *(synonyms.get("tissues", []) or []),
            ]
        )
        free_terms = self._unique(extracted.get("free_terms", []) or [])

        return ClassificationContext(
            organisms=organisms,
            strategies=strategies,
            conditions=conditions,
            tissues=tissues,
            free_terms=free_terms,
        )

    def _classify_with_ollama(
        self,
        source: str,
        record: Dict[str, Any],
        context: ClassificationContext,
    ) -> Optional[Dict[str, Any]]:
        if not self.ollama_client:
            return None

        prompt = self._build_ollama_prompt(source=source, record=record, context=context)

        try:
            raw = self.ollama_client.generate(prompt, temperature=0.0)
            match = re.search(r"\{.*\}", raw, re.DOTALL)
            if not match:
                return None
            parsed = json.loads(match.group(0))
            return self._validate_classification(parsed, model_source="ollama")
        except Exception:
            return None

    def _build_ollama_prompt(self, source: str, record: Dict[str, Any], context: ClassificationContext) -> str:
        compact_record = {
            "source": source,
            "title": record.get("title", ""),
            "abstract": record.get("abstract", ""),
            "organism": record.get("organism", ""),
            "project_type": record.get("project_type", ""),
            "publication_type": record.get("publication_type", ""),
            "description": record.get("description", ""),
        }

        return (
            "Eres un clasificador cientifico de resultados de busqueda. "
            "Responde SOLO JSON valido con estas llaves: "
            "relevance_label, relevance_score, reason_short, tags, evidence_level, model_source.\n\n"
            "Reglas:\n"
            "- relevance_label: alta, media, baja\n"
            "- relevance_score: numero entre 0 y 1\n"
            "- reason_short: 1 o 2 frases en espanol\n"
            "- tags: lista de strings\n"
            "- evidence_level: directa, indirecta, débil\n"
            "- model_source: ollama\n\n"
            f"Contexto consulta: {json.dumps(context.__dict__, ensure_ascii=False)}\n"
            f"Resultado: {json.dumps(compact_record, ensure_ascii=False)}"
        )

    def _heuristic_pubmed(self, record: Dict[str, Any], context: ClassificationContext) -> Dict[str, Any]:
        text = " ".join(
            [
                str(record.get("title", "")),
                str(record.get("abstract", "")),
                str(record.get("journal", "")),
                str(record.get("publication_type", "")),
            ]
        )

        score, tags, evidence_level = self._score_text(text=text, record=record, context=context, source="pubmed")
        label = self._label_from_score(score)

        if label == "alta":
            reason = "Alta concordancia con organismo/condicion/estrategia de la consulta."
        elif label == "media":
            reason = "Coincidencia parcial con la consulta; revisar detalle biologico."
        else:
            reason = "Baja concordancia con los criterios principales de la consulta."

        return {
            "relevance_label": label,
            "relevance_score": score,
            "reason_short": reason,
            "tags": tags,
            "evidence_level": evidence_level,
            "model_source": "heuristic",
        }

    def _heuristic_bioproject(self, record: Dict[str, Any], context: ClassificationContext) -> Dict[str, Any]:
        hierarchy = record.get("sra_hierarchy", {}) or {}
        strategy_tokens = self._collect_sra_strategies(hierarchy)

        text = " ".join(
            [
                str(record.get("title", "")),
                str(record.get("description", "")),
                str(record.get("organism", "")),
                " ".join(strategy_tokens),
            ]
        )

        score, tags, evidence_level = self._score_text(text=text, record=record, context=context, source="bioproject")
        label = self._label_from_score(score)

        if strategy_tokens:
            tags.extend([f"estrategia:{s}" for s in strategy_tokens[:2] if f"estrategia:{s}" not in tags])

        if label == "alta":
            reason = "Proyecto fuertemente alineado con la consulta y metadatos SRA disponibles."
        elif label == "media":
            reason = "Proyecto potencialmente util; requiere revisar experimentos y muestras."
        else:
            reason = "Proyecto con baja alineacion frente a la pregunta planteada."

        return {
            "relevance_label": label,
            "relevance_score": score,
            "reason_short": reason,
            "tags": tags,
            "evidence_level": evidence_level,
            "model_source": "heuristic",
        }

    def _score_text(
        self,
        text: str,
        record: Dict[str, Any],
        context: ClassificationContext,
        source: str,
    ) -> Tuple[float, List[str], str]:
        norm_text = self._normalize(text)
        tags: List[str] = []

        org_matches = self._find_matches(norm_text, context.organisms)
        cond_matches = self._find_matches(norm_text, context.conditions)
        strat_matches = self._find_matches(norm_text, context.strategies)
        tissue_matches = self._find_matches(norm_text, context.tissues)
        free_matches = self._find_matches(norm_text, context.free_terms)

        score = 0.0
        score += 0.35 * self._coverage(org_matches, context.organisms)
        score += 0.20 * self._coverage(strat_matches, context.strategies)
        score += 0.20 * self._coverage(cond_matches, context.conditions)
        score += 0.10 * self._coverage(tissue_matches, context.tissues)
        score += 0.10 * self._coverage(free_matches, context.free_terms)

        completeness = self._field_completeness(record=record, source=source)
        score += 0.05 * completeness

        score = round(min(score, 1.0), 2)

        tags.extend([f"organismo:{m}" for m in org_matches[:2]])
        tags.extend([f"condicion:{m}" for m in cond_matches[:2]])
        tags.extend([f"estrategia:{m}" for m in strat_matches[:2]])
        tags.extend([f"tejido:{m}" for m in tissue_matches[:2]])

        publication_type = self._normalize(str(record.get("publication_type", "")))
        if source == "pubmed" and "review" in publication_type:
            tags.append("evidencia:review")

        if not tags:
            tags.append("alineacion:baja")

        label = self._label_from_score(score)
        evidence_level = "indirecta"
        if source == "pubmed" and "review" in publication_type:
            evidence_level = "indirecta" if label != "baja" else "débil"
        elif label == "alta" and (org_matches or strat_matches):
            evidence_level = "directa"
        elif label == "media":
            evidence_level = "indirecta"
        else:
            evidence_level = "débil"

        return score, self._unique(tags), evidence_level

    def _collect_sra_strategies(self, hierarchy: Dict[str, Any]) -> List[str]:
        strategies: List[str] = []
        for sample in hierarchy.values():
            experiments = sample.get("experiments", []) if isinstance(sample, dict) else []
            for exp in experiments:
                metadata = exp.get("metadata", {}) if isinstance(exp, dict) else {}
                strategy = metadata.get("library_strategy")
                if strategy:
                    strategies.append(str(strategy))
        return self._unique(strategies)

    def _field_completeness(self, record: Dict[str, Any], source: str) -> float:
        if source == "pubmed":
            fields = ["pmid", "title", "abstract", "journal", "year"]
        else:
            fields = ["bioproject", "title", "description", "organism", "sra_experiments_count"]

        valid = 0
        for field in fields:
            value = record.get(field)
            if value is None:
                continue
            if isinstance(value, str) and value.strip().upper() in {"", "NA"}:
                continue
            valid += 1

        return valid / max(len(fields), 1)

    def _validate_classification(self, payload: Dict[str, Any], model_source: str) -> Optional[Dict[str, Any]]:
        try:
            label = str(payload.get("relevance_label", "")).strip().lower()
            if label not in VALID_LABELS:
                return None

            score = float(payload.get("relevance_score", 0))
            score = max(0.0, min(1.0, score))

            reason = str(payload.get("reason_short", "")).strip()
            if not reason:
                return None

            evidence = str(payload.get("evidence_level", "")).strip().lower()
            if evidence == "debil":
                evidence = "débil"
            if evidence not in VALID_EVIDENCE:
                return None

            tags = payload.get("tags", [])
            if not isinstance(tags, list):
                tags = []
            tags = [str(tag).strip() for tag in tags if str(tag).strip()]

            return {
                "relevance_label": label,
                "relevance_score": round(score, 2),
                "reason_short": reason,
                "tags": self._unique(tags) or ["alineacion:modelo"],
                "evidence_level": evidence,
                "model_source": model_source,
            }
        except Exception:
            return None

    @staticmethod
    def _label_from_score(score: float) -> str:
        if score >= 0.75:
            return "alta"
        if score >= 0.45:
            return "media"
        return "baja"

    @staticmethod
    def _coverage(matches: List[str], terms: List[str]) -> float:
        if not terms:
            return 0.0
        return min(len(matches) / len(terms), 1.0)

    @staticmethod
    def _normalize(value: str) -> str:
        value = unicodedata.normalize("NFD", value.lower())
        value = "".join(ch for ch in value if unicodedata.category(ch) != "Mn")
        return re.sub(r"\s+", " ", value).strip()

    def _find_matches(self, norm_text: str, terms: List[str]) -> List[str]:
        matches: List[str] = []
        for term in terms:
            token = self._normalize(str(term))
            if token and token in norm_text:
                matches.append(str(term))
        return self._unique(matches)

    @staticmethod
    def _unique(values: List[Any]) -> List[str]:
        seen = set()
        ordered: List[str] = []
        for value in values:
            if value is None:
                continue
            text = str(value).strip()
            if not text:
                continue
            key = text.lower()
            if key in seen:
                continue
            seen.add(key)
            ordered.append(text)
        return ordered

# BioProject Metadata: Current vs Potential

## 📊 INFORMACIÓN ACTUAL QUE ESTAMOS OBTENIENDO

Actualmente el fetcher extrae **4 campos** de cada BioProject:

| Campo | Ejemplo | Utilidad |
|-------|---------|----------|
| **sra_id** | 1220320 | Identificador SRA vinculado al proyecto |
| **bioproject** | PRJNA1220320 | Identificador único del BioProject en NCBI |
| **title** | "Nicotiana benthamiana Raw sequence reads" | Título descriptivo del proyecto |
| **description** | "This study aimed at exploring..." | Descripción detallada del proyecto |
| **organism** | (vacío) | Organismo estudiado |
| **fetched_at** | 2026-02-12T11:27:18 | Timestamp de la búsqueda |

---

## 🔍 INFORMACIÓN DISPONIBLE EN LA API DEL BIOPROJECT QUE NO ESTAMOS USANDO

### Metadatos de Organismo
- **Project_Organism_Name**: Nombre científico completo del organismo
- **Project_Organism_Strain**: Cepa específica del organismo
- **Project_Organism_Breed**: Raza/variante del organismo
- **Project_Organism_Cultivar**: Cultivar de plantas

**Utilidad**: Permitiría filtrar búsquedas por cepa específica, cultivares, etc.

### Información del Proyecto
- **Project_Type**: Tipo de proyecto (umbrella, single, organization)
- **Project_Method**: Metodología del proyecto
- **Project_Category**: Categoría del proyecto
- **Relevance**: Evaluación de relevancia del proyecto

**Utilidad**: Clasificar proyectos por tipo y metodología

### Información de Diseño Experimental
- **Experimental_Design**: Descripción del diseño experimental
- **Publications**: Publicaciones asociadas
- **Locus_Tag_Prefix**: Prefijo de etiquetas de locus (para anotaciones)

**Utilidad**: Entender rápidamente el enfoque experimental sin leer toda la descripción

### Datos Temporales
- **Submission_Date**: Fecha original de envío del proyecto
- **Public_Date**: Fecha cuando se hizo público
- **Update_Date**: Última fecha de actualización
- **Project_Grant_Title**: Título de la beca que financió el proyecto
- **Project_Grant_Agency**: Agencia que financió (NIH, NSF, etc.)

**Utilidad**: Conocer antigüedad del proyecto y financiamiento

### Información de Acceso y Datos
- **Total_Studies**: Número de estudios en el proyecto
- **Total_Experiments**: Número de experimentos
- **Total_Runs**: Número de secuenciaciones/runs
- **Total_Bases**: Número total de bases secuenciadas
- **Database_Records**: Cantidad de registros en diferentes bases de datos
- **Data_Size**: Tamaño total de datos (en GB)

**Utilidad**: Evaluar escala del proyecto, cantidad de datos disponibles

### Información de Mantenimiento
- **Project_Owner_Mail**: Email del propietario del proyecto
- **Project_Owner_Organization**: Institución del investigador principal

**Utilidad**: Contacto directo con investigadores

### Información de Temas/Keywords
- **Project_Keyword**: Palabras clave del proyecto
- **Disease**: Enfermedades asociadas al estudio
- **Environmental_Sample**: Si es muestra ambiental
- **Isolate**: Si es un aislado

**Utilidad**: Clasificación automática, búsquedas más específicas

---

## 💡 RECOMENDACIONES POR PRIORIDAD

### 🔴 ALTA PRIORIDAD (Muy útil, bajo esfuerzo)
1. **Project_Organism** - Ya existe en la respuesta, solo hay que extraerlo
2. **Submission_Date, Public_Date** - Para tracking temporal
3. **Total_Studies, Total_Experiments, Total_Runs** - Para evaluar escala del proyecto
4. **Project_Grant_Agency** - Para identificar financiamiento

**Impacto**: Permitiría filtrar por organism, fecha, y tamaño del proyecto

### 🟡 MEDIA PRIORIDAD (Útil, moderado esfuerzo)
1. **Project_Organism_Strain** - Para searches específicas por cepa
2. **Publications** - Para validar relevancia
3. **Experimental_Design** - Más información sin leer descripción completa
4. **Data_Size, Total_Bases** - Para evaluar viabilidad de análisis

**Impacto**: Mejor clasificación y evaluación de proyectos

### 🟢 BAJA PRIORIDAD (Niche, puede dejarse para después)
1. **Project_Owner_Mail** - Para contacto directo (privacidad)
2. **Project_Type** - Clasificación general
3. **Locus_Tag_Prefix** - Principalmente para anotaciones de genomas

**Impacto**: Información adicional para casos específicos

---

## 📝 SUGERENCIA DE NUEVA ESTRUCTURA CSV

```
sra_id,bioproject,title,description,organism,strain,submission_date,public_date,
total_studies,total_runs,total_bases,data_size_gb,grant_agency,fetched_at
```

### Ejemplo de fila mejorada:
```
1220320,PRJNA1220320,"Nicotiana benthamiana Raw sequence reads",...,
Nicotiana benthamiana,Col-0,2023-05-12,2023-06-01,
1,24,2.8e9,15.2,NSF,2026-02-12
```

---

## 🔧 IMPLEMENTACIÓN SUGERIDA

1. **Fase 1** (INMEDIATA): Agregar organism, strain, grant_agency, total_runs
2. **Fase 2** (PRÓXIMA): Agregar dates (submission, public), total_bases, data_size
3. **Fase 3** (FUTURA): Agregar experimental_design, publications, owner_contact

Esto permitiría:
- ✅ Filtrar por organismo y cepa específica
- ✅ Evaluar rápidamente escala de proyectos (runs, bases)
- ✅ Conocer fuentes de financiamiento
- ✅ Rastrear evolución temporal de proyectos
- ✅ Validar proyectos por publicaciones asociadas

#!/usr/bin/env python3
"""
Quick test to verify metadata extraction works correctly
"""

# Simulate NCBI BioProject API response
sample_bioproject = {
    "Project_Acc": "PRJNA1234567",
    "Project_Title": "Test Project",
    "Project_Description": "Test description",
    "Project_Organism": "Arabidopsis thaliana",
    "Project_OrganismList": [
        {
            "OrganismName": "Arabidopsis thaliana",
            "Strain": "Col-0",
            "Cultivar": "Columbia"
        }
    ],
    "Submission_Date": "2024-01-15",
    "Public_Date": "2024-02-01",
    "Total_Studies": "5",
    "Total_Runs": "50",
    "Publications": [
        {
            "PMID": "12345678",
            "Title": "Test Publication 1"
        },
        {
            "PMID": "87654321",
            "Title": "Test Publication 2"
        }
    ],
    "Experimental_Design": "Study of drought response in plants"
}

# Simulate the extraction logic
def extract_bioproject_metadata(summary):
    """Extract metadata from BioProject summary response."""
    try:
        bioproject = summary.get("Project_Acc", "")
        title = summary.get("Project_Title", "")
        description = summary.get("Project_Description", "")
        
        # Organismos
        organism = summary.get("Project_Organism", "")
        organism_list = summary.get("Project_OrganismList", [])
        organisms = organism
        
        if organism_list and isinstance(organism_list, (list, dict)):
            if isinstance(organism_list, list):
                organisms_found = []
                for org in organism_list:
                    if isinstance(org, dict) and org.get("OrganismName"):
                        organisms_found.append(org.get("OrganismName", ""))
                if organisms_found:
                    organisms = "; ".join(organisms_found)
        
        # Strain y Cultivar
        strain = ""
        cultivar = ""
        if organism_list and isinstance(organism_list, list):
            for org in organism_list:
                if isinstance(org, dict):
                    if org.get("Strain") and not strain:
                        strain = org.get("Strain", "")
                    if org.get("Cultivar") and not cultivar:
                        cultivar = org.get("Cultivar", "")
        
        # Fechas
        submission_date = summary.get("Submission_Date", "")
        public_date = summary.get("Public_Date", "")
        
        # Conteos
        total_studies = summary.get("Total_Studies", "")
        total_runs = summary.get("Total_Runs", "")
        
        # Publicaciones
        publications_list = summary.get("Publications", [])
        publications = ""
        if publications_list and isinstance(publications_list, list):
            publications_found = []
            for pub in publications_list:
                if isinstance(pub, dict):
                    pmid = pub.get("PMID", "")
                    if pmid:
                        publications_found.append(f"PMID:{pmid}")
            if publications_found:
                publications = "; ".join(publications_found)
        
        # Diseño experimental
        experimental_design = summary.get("Experimental_Design", "")
        
        if not bioproject:
            return None
        
        return {
            "bioproject": bioproject,
            "title": title,
            "description": description,
            "organisms": organisms,
            "strain": strain,
            "cultivar": cultivar,
            "submission_date": submission_date,
            "public_date": public_date,
            "total_studies": total_studies,
            "total_runs": total_runs,
            "publications": publications,
            "experimental_design": experimental_design
        }
    except Exception as e:
        print(f"Error extracting metadata: {e}")
        return None

# Test the extraction
print("=" * 70)
print("Testing BioProject Metadata Extraction")
print("=" * 70)

result = extract_bioproject_metadata(sample_bioproject)

if result:
    print("\n✅ Extraction successful!")
    print("\nExtracted fields:")
    for key, value in result.items():
        if value:
            print(f"  {key:20s}: {value}")
        else:
            print(f"  {key:20s}: (empty)")
else:
    print("\n❌ Extraction failed!")

print("\n" + "=" * 70)
print("Expected CSV header:")
print("sra_id,bioproject,title,organisms,strain,cultivar,submission_date,public_date,total_studies,total_runs,publications,experimental_design,description,fetched_at")
print("=" * 70)

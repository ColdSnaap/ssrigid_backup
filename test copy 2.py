from mp_api.client import MPRester
with MPRester(api_key="DAtVcNKmI2Bg2smQxHWh81PQlJEhcL13") as mpr:
    data = mpr.materials.search(material_ids=["mp-2646995"])

material_ids = ["mp-2646995"]
for doc in data:
    structure = doc.structure
    material_id = doc.material_id

    structure.to("POSCAR", "posacr_test")
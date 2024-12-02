import pathlib
from astropy.table import Table, vstack
from astroquery.simbad import Simbad

base_path = pathlib.Path("simbad_query_results")
base_path.mkdir(parents=True, exist_ok=True)
max_hip_id = 120416
query_batch_size = 2000

custom_simbad = Simbad()
custom_simbad.TIMEOUT = 99999
custom_simbad.add_votable_fields(
    "ids",
    "ra(d2000)",
    "dec(d2000)",
    "pmra",
    "pmdec",
    "plx",
    "plx_error",
    "flux(U)",
    "flux(B)",
    "flux(V)",
    "flux(R)",
    "flux(I)",
    "flux(G)",
    "flux(J)",
    "flux(H)",
    "flux(K)",
    "otype(N)",
    "sptype",
    "rv_value",
    "rvz_error",
)

# dont make a giant query, make few smaller queries
for batch in range(1, max_hip_id // query_batch_size + 1):
    hip_id = ["HIP " + str(i) for i in range(1 + (batch - 1) * query_batch_size, 1 + batch * query_batch_size)]
    result = custom_simbad.query_objects(hip_id)
    result.write(base_path / f"simbad_hipparcos_{str(batch)}.dat", format="ascii", overwrite=True)

# do the remaining id
hip_id = ["HIP " + str(i) for i in range(batch * query_batch_size, max_hip_id + 1)]
result = custom_simbad.query_objects(hip_id)
result.write(base_path / f"simbad_hipparcos_{str(batch + 1)}.dat", format="ascii", overwrite=True)

# merge all the tables
files_list = base_path.glob("simbad_*.dat")
table_list = []

for file in files_list:
    table_list.append(Table.read(file, format="ascii"))

simbad_table = vstack(table_list)
simbad_table.write(base_path / "hip_simbad.dat", format="ascii")

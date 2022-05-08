import subprocess, os
from fireworks import (
    FiretaskBase,
    FWAction,
    explicit_serialize,
    Workflow,
    Firework,
    LaunchPad,
)


@explicit_serialize
class OptimizeComplex(FiretaskBase):
    required_params = ["input_path"]
    optional_params = []

    def run_task(self, fw_spec):

        # Variables
        orca_inp_file = self["input_path"]
        print("It is running this codebase!")
        run_orca = f"/opt/orca-5.0.2/orca {orca_inp_file} '--map-by hwthread' > {orca_inp_file.split('/')[-1].split('.')[0]}_sp_orca.out"
        status_code = subprocess.call(run_orca, stdout=subprocess.PIPE, shell=True)
        return FWAction(stored_data={"output": "some energy", "status_code": status_code})


if __name__ == "__main__":
    #inputs = [
    #    "/home/jovyan/shared_scratch/individual_xyz_data/inputs/LEKRID.xyz.inp",
    #    "/home/jovyan/shared_scratch/individual_xyz_data/inputs/OXACUN.xyz.inp",
    #]
#    inputs = ["/home/jovyan/oxo_intermediates/inputs/ERAPUI.xyz.inp",
#            "/home/jovyan/oxo_intermediates/inputs/FEKGUA.xyz.inp"
#            ]
    molecules = open('oxo_int_charge_mult_success.txt','r').readlines()
    launchpad = LaunchPad(
        host="localhost",
        name="fw_oal",
        port=27017,
        username="fw_oal_admin",
        password="gfde223223222rft3",
    )
    fws = []
    for mol in molecules[:900]:
        mol_name = mol.split('\n')[0]
        input_ = f"/home/jovyan/oxo_intermediates/inputs/{mol_name}.xyz_charge_mult.inp"
        fw = Firework(
            OptimizeComplex(input_path=f"{input_}"), name=f"lean_orca_FW_{input_}"
        )
        fws.append(fw)
    wf = Workflow(fws, name="lean_orca_func")
    launchpad.add_wf(wf)
    print("Executing")

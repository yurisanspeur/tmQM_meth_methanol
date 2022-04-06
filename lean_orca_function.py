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
        process = subprocess.Popen(run_orca.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        return FWAction(stored_data={"output": "some energy"})


if __name__ == "__main__":
    inputs = [
        "/home/jovyan/shared_scratch/individual_xyz_data/inputs/LEKRID.xyz.inp",
        "/home/jovyan/shared_scratch/individual_xyz_data/inputs/OXACUN.xyz.inp",
    ]
    launchpad = LaunchPad(
        host="localhost",
        name="fw_oal",
        port=27017,
        username="fw_oal_admin",
        password="gfde223223222rft3",
    )
    fws = []
    for input_ in inputs:
        fw = Firework(
            OptimizeComplex(input_path=f"{input_}"), name=f"lean_orca_FW_{input_}"
        )
        fws.append(fw)
    wf = Workflow(fws, name="lean_orca_func")
    launchpad.add_wf(wf)
    print("Executing")

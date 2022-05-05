find . -name "*.xyz" -printf "%P\n" | while read -r test; do
  echo $test
  var=$(echo $test | awk -F '.' '{print $1".xyz"}')
  echo $var
  match=$(echo $test | awk -F '.' '{print $1}') # without the extension
  echo $match
#  echo "s/ABAFUF/$test/g"
#  cat /home/jovyan/shared_scratch/individual_xyz_data/run_orca_parallel_test/b3lyp/pbe0/OXAFID_sp_orca.inp | sed "s/AAZDCO.xyz/${var}_oxo_intermediate.xyz/g" > "inputs/${var}.inp"
  #cat "$test" | grep MND | awk '{print $3","$5","$6","$7","$9","$10","$11","$17","$18","$19}'
  # Get the mult and charge from the csv file obtained from python script
  mult=$(cat ./charge_mult_oxo_success.csv | grep $match | awk -F ',' '{print $NF}')
  echo $mult
  # Mutate the charge and the multiplicity based on the intermediate
  if [ ! -z $mult ]; then
    sed "s/\s1/ ${mult}/g" "inputs/${var}.inp" > "inputs/${var}_mult.inp"
  fi
#  # Hardcoded charge of -1 for the oxo_intermediate
  sed "s/xyzfile\s0/xyzfile -1/g" "inputs/${var}_mult.inp" > "inputs/${var}_charge_mult.inp"
done








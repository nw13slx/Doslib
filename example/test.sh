mkdir run 
cd run
../../doslib/vasp_dos -path ../inputfile/ -p 0 -f 1 --peratom --write_pdos 
cd ../

# Individual-project
Protocol for reproducing the report results

## 1. Environment building
### 1) AlphaFold3: 
build the environment following the AlphaFold3 GitHub (https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md)

### 2) MemProtMD for membrane insertion: 
using Google Colab website (https://colab.research.google.com/github/pstansfeld/MemProtMD/blob/main/MemProtMD_Insane.ipynb#scrollTo=0GnUhNBFAa3t) or add ... to .bashrc to build the environment called memprotmd.

## 2. Structure prediction
   
### 1) Sequence obtained: 
protein sequence (E.coli) obtained from UniProt (https://www.uniprot.org/); small molecule will get the SMILES string from chEBI (https://www.ebi.ac.uk/chebi/init.do) or Protein Data Bank (PDB) (https://www.rcsb.org/).

### 2) Predicting structures

#### AlphaFold3

Create a JSON input file (`input.json`), then run:

**Command**
```
python run_alphafold.py \
  --json_path=/root/af_input/input.json \
  --model_dir=/root/models \
  --output_dir=/root/af_output
```
Example `input.json`

```
{
  "name": "WecA-glc",
  "sequences": [
    {
      "protein": {
        "id": ["A"],
        "sequence": "MNLLTVSTDLISIFLFTTLFLFFARKVAKKVGLVDKPNFRKRHQGLIPLVGGISVYAGICFTFGIVDYYIPHASLYLACAGVLVFIGALDDRFDISVKIRATIQAAVGIVMMVFGKLYLSSLGYIFGSWEMVLGPFGYFLTLFAVWAAINAFNMVDGIDGLLGGLSCVSFAAIGMILWFDGQTSLAIWCFAMIAAILPYIMLNLGILGRRYKVFMGDAGSTLIGFTVIWILLETTQGKTHPISPVTALWIIAIPLMDMVAIMYRRLRKGMSPFSPDRQHIHHLIMRAGFTSRQAFVLITLAAALLASIGVLAEYSHFVPEWVMLVLFLLAFFLYGYCIKRAWKVARFIKRVKRRLRRNRGGSPNLTK"
      }
    },
    {
      "ligand": {
        "id": ["L"],
        "smiles": "	CC(=O)N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]1OP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1ccc(=O)[nH]c1=O"
      }
    }
  ],
  "modelSeeds": [1],
  "dialect": "alphafold3",
  "version": 1
}
```

#### Chai-1: put the same sequences in the Chai-1 website (https://www.chaidiscovery.com/).

## 3. Generate topology and energy minimization before inserting the membrane:

  ```
  gmx pdb2gmx -f weca-glc.pdb -o weca-glc-p2g.pdb
  gmx editconf -f weca-glc-p2g.pdb -o weca-glc-box.pdb -d 2 -c
  gmx grompp -f em100.mdp -c weca-glc-box.pdb -o weca-glc-em -maxwarn 2
  gmx mdrun -deffnm weca-glc-em -v -ntmpi 1 -c weca-glc-em.pdb
  ```
  
## 4. Membrane insertion

   If you use Google Colab for MemprotMD, upload the protein sequence, change the membrane type to 'POPE 7 POPG 2 CARD 1', and then run all. If you use the MemprotMD locally, run the following command after activating the memprotmd environment:

   ```
   module load GCC/11.2.0 Boost/1.77
   python pdb-fix.py weca-glc-em.pdb 
   python memprotmd-insane-martini3.py fixed-weca-glc-em.pdb '-l POPE:7 -l POPG:2 -l CARD:1'
   python cg2at -a weca-glc-em.pdb -c equil2-1.gro -loc CG2AT -w tip3p -fg martini_3-0_charmm36
   ```

   Reminder: the MemprotMD can only insert the membrane into the individual protein, not the protein-ligand complex. Therefore, you need to align the structure after insertion and the protein-ligand complex (e.g., `weca-glc-em.pdb`).

   **Align the structure in Pymol**

   ```
   pymol weca-glc-em.pdb final_cg2at_aligned.pdb
 
   In pymol:
	 align weca-glc-em, final_cg2at_aligned and pol
	 save weca-glc-aligned.pdb, weca-glc-em
	 remove pol
	 save AT-system.pdb, final_cg2at_aligned
   ```

## 5. 10ns NPT equilibration under positional restraints

## 6. Production MD run


**The MD run will cost more than hours so that the trajectories will be provided in ....**

## 7. Data analysis and post-processing
   

   



   

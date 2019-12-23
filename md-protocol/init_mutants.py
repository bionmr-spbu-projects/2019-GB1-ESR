#!/usr/bin/python2.7

#
# WARNING: No other global imports are allowed here, use imports inside class functions instead
#          Current implementation can not track such dependencies and your run will fail

from runmd2.MD import MD


class MD_1pgb(MD):

    def __init__(self, init_struct ):
        import os
        from pyxmol import readPdb
        restricted_structure = readPdb(init_struct)[0]
        wd = os.path.splitext(os.path.basename(init_struct))[0]
        super(MD_1pgb, self).__init__(name=wd, trj_home_dir=wd)

        restricted_structure.writeAsPdb(open("input.pdb", "w"))
        self.save_state("saved_state.pickle")
        self.save_state("auto_dump.pickle")

    def run_setup(self):
        import os
        from pyxmol import readPdb
        from pyxmol.predicate import rId

        self.log("Setup > ")
        self.tleaprc.source(os.environ["AMBERHOME14"] + "/dat/leap/cmd/leaprc.ff14SB")
        self.tleaprc.add_command("loadamberparams frcmod.ionsjc_tip3p")
        self.tleaprc.add_command("loadoff /home/seva/GB1/from_purcell/cml.lib")
        self.tleaprc.add_command("loadoff /home/seva/GB1/from_purcell/ntcml.lib")
        self.tleaprc.add_command("loadamberparams /home/seva/GB1/from_purcell/cml.frcmod")
        self.tleaprc.load_pdb("input.pdb")
        self.tleaprc.solvate_oct("TIP3PBOX", 10.0)
        # self.tleaprc.add_command("solvatebox wbox TIP3PBOX 10")
        self.tleaprc.add_ions("Na+",target_charge=0)
        # self.tleaprc.add_command("addIonsRand wbox Na+ 0")
        self.tleaprc.save_params(output_name=MD._build_dir+"/box")
        self.tleaprc.save_pdb(output_name=MD._build_dir+"/box")
        self.tleaprc.add_command("quit")
        
        self.min1_parameters.set(
            imin=1,
            maxcyc=500,
            ncyc=100,
            ntb=1,
            ntr=1,
            cut=10.0
        ).add_atom_pin(200, None, [(1,56)])
        self.min2_parameters.set(
            imin=1,
            maxcyc=100,
            ncyc=200,
            ntb=1,
            ntr=0,
            cut=10.0
        ).add_atom_pin(200, None, [(1,56)])
        self.heat_parameters.set(
            imin=0,
            irest=0,
            ntx=1,
            ntb=1,
            ntr=1,
            ntc=2,
            ntf=2,
            tempi=0.0,
            temp0=293.0,
            ntt=1,
            nstlim=10000,
            dt=0.001,
            ntpr=50,
            ntwx=50,
            ntwr=50,
            ioutfm=1
        ).add_atom_pin(10, None, [(1,56)])
        self.equil_parameters.set(
            imin=0,
            irest=1,
            ntx=5,
            ntb=2,
            iwrap=1,
            ntt=3,
            gamma_ln=2.0,
            ig=-1,
            tempi=293.0,
            temp0=293.0,
            ntp=1,
            pres0=1.0,
            taup=2.0,
            cut=10.0,
            ntr=1,
            ntc=2,
            ntf=2,
            nstlim=500000,
            dt=0.002,
            ntpr=500,
            ntwx=500,
            ntwr=500000,
            ioutfm=1
        )
        self.run_parameters.set(
            imin=0,
            irest=1,
            ntx=5,
            ntb=2,
            iwrap=1,
            ntt=3,
            gamma_ln=2.0,
            tempi=293.0,
            temp0=293.0,
            ntp=1,
            pres0=1.0,
            taup=2.0,
            cut=10.5,
            ntr=0,
            ntc=2,
            ntf=2,
            nstlim=500000,
            dt=0.002,
            ntpr=500,
            ntwx=500,
            ntwr=500000,
            ioutfm=1
        )

        self._pmemd_executable = ["/home/sergei/amber16/bin/pmemd.cuda"]

        self.build()

        self.restricted_structure = readPdb(self.tleaprc.pdb_output_name+".pdb")[0].asResidues >> (rId<=56)

        self.minimize()
        self.heat()
        self.equilibrate()

        self.setup_is_done = True
        # self._pmemd_executable = ["mpirun","-np","2","/home/sergei/amber16/bin/pmemd.cuda.MPI"]
        self.save_state("saved_state.pickle")
        self.log("Setup < ")

    def run_continue(self):
        import random


        self.log("Continue > ")

        self.log("")
        self.log("Unconstrained simulation")
        self.log("")
        
        ats = self.restricted_structure.asAtoms
        
        while self.current_step < 20000:
            self.run_parameters["ig"] = random.randint(1,100000)
            self.log("Random seed is %d " % self.run_parameters["ig"])
            self.do_md_step()
            self.put_frame(ats, -1, -1)
            ats.writeAsPdb(open(MD._run_dir+"/run%05d.pdb" % self.current_step, "w"))
        
        self.save_state("all_done.pickle")

        self.log("Continue < ")


import os
import glob
wd = os.path.abspath(".")
for struct in glob.glob("initial_structures/????_*C_CML.pdb"):
    os.chdir(wd)
    MD_1pgb(struct)




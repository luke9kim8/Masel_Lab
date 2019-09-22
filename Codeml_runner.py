from Bio.Phylo.PAML import codeml
import subprocess

class Codeml_runner:
    def __init__(self):
        self
    def run_bash_command(self,cmd):
        """
        Constructor for Codeml_runner obj
        :param cmd:
        :return:
        """
        subprocess.Popen(cmd, shell=True, executable='/bin/bash')

    def run_codeml(self, seqfile, treefile):
        """
        Writes a domain tree with evolution times; The tree is stored in scratch directory
        :param seqfile:
        :param outfile:
        :param treefile:
        :return:
        """
        cml = codeml.Codeml()
        cml.read_ctl_file("Codeml/codeml.ctl")
        cml.tree = treefile
        cml.out_file = 'Codeml/Codeml_files/codeml.file'
        cml.alignment = seqfile
        cml.write_ctl_file()
        self.run_bash_command('Codeml/codeml codeml.ctl')

#!/usr/bin/env python3
import os
import logging
import subprocess as sp
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
from Bio.Seq import Seq


def read_gene_families(gftxt, protfile = None, cdsfile = None, wrkdir = None):
    """
    Read gene families and sequences and put them into the GeneFamily Class
    """
    gene_families = []
    if protfile is None and cdsfile is None:
        logging.info("Gene families need to have sequences!")
        with open(gftxt, 'r') as f:
            for line in f:
                line = line.rstrip()
                x = line.split()
                gf_id = x.pop(0)[:-1]
                gf_genes = x
                gene_families.append(GeneFamily(gf_id=gf_id, gf_members=gf_genes))
        return gene_families
    
    if protfile is not None:
        prot = SeqIO.to_dict(SeqIO.parse(protfile, "fasta"))

    if cdsfile is not None:
        cds = SeqIO.to_dict(SeqIO.parse(cdsfile, "fasta"))
    
    with open(gftxt, 'r') as handle:
        for line in handle:
            line = line.rstrip()
            x = line.split()
            gf_id = x.pop(0)[:-1]
            gf_genes = x
            gf_prot = {}
            gf_cds = {}
            for gid in x:
                if prot[gid][-1:].seq == '*':
                    gf_prot[gid] = prot[gid][:-1]
                else:
                    gf_prot[gid] =prot[gid]
                if cds[gid][-3:].seq == "TAA" or \
                cds[gid][-3:].seq == "TAG" or \
                cds[gid][-3:].seq == "TGA":
                    gf_cds[gid] = cds[gid][:-3]
                else:
                    gf_cds[gid] = cds[gid]
            gene_families.append(GeneFamily(gf_id = gf_id, gf_members = gf_genes, 
                prot_seqs = gf_prot, cds_seqs = gf_cds, wrkdir=wrkdir))
    return gene_families

def split_gene_families(gftxt, max=200, min=4, num=2):
    gfs = read_gene_families(gftxt)
    gfsets = [[None] * 1 for i in range(num)] 
    for i in range(len(gfs)):
        l = len(gfs[i].members)
        if l < min or l > max:
            continue
        else:
            m = i % num
            fname = 'GeneFamilySet' + str(m) + '.txt'
            with open(fname, "a") as out:
                out.write('%s: %s\n' % (gfs[i].id, ' '.join(gfs[i].members)))
    
def _write_fasta(fname, seq_dict):
    with open(fname, "w") as f:
        for i, s in seq_dict.items():
            #print(i)
            #print(s.seq)
            f.write(">%s\n%s\n" % (i, s.seq))
    return fname

def _mkdir(dirname):
    if os.path.isdir(dirname):
        logging.warning("dir {} exists!".format(dirname))
    else:
        os.makedirs(dirname)
    return dirname
    
class GeneFamily(object):
    def __init__(self, gf_id, gf_members, prot_seqs=None, cds_seqs=None,
            #prot_aln = None, cds_aln = None,
            wrkdir = None, phylo = None, 
            aligner = "muscle", phylotool = "raxml"):
        self.id = gf_id
        self.members = gf_members
        self.prot = prot_seqs
        self.cds = cds_seqs
        #self.prot_aln = prot_aln
        #self.cds_aln = cds_aln
        self.aligner = aligner
        self.phylotool = phylotool
        if wrkdir == None:
            self.wrkdir = os.path.join('run', gf_id)
        else:
            self.wrkdir = os.path.join(wrkdir, gf_id)
    
    def get_wrkdir(self):
        print(self.wrkdir)

    def get_id(self):
        print(self.id)
    
    def get_prot(self):
        #print(self.id)
        pep = self.id + ".pep"
        for m in self.members:
            #print(m)
            #print(self.prot[m].seq)
            print(self.prot[m].format("fasta").rstrip())

    def get_cds(self):
        #print(self.id)
        cds = self.id + ".cds"
        for m in self.members:
            #print(m)
            #print(self.cds[m].seq)
            print(self.cds[m].format("fasta").rstrip())
    
    def obtain_msa(self):
        wd = _mkdir(self.wrkdir)
        fname = os.path.join(self.wrkdir, self.id + ".pep")
        _write_fasta(fname, self.prot)
        fname = os.path.join(wd, self.id + ".cds")
        _write_fasta(fname, self.cds)
        if (self.aligner == "muscle"):
            self.run_muscle()
        elif (self.aligner == "mafft"):
            pass
            #run_mafft()
        self.run_trimal()

    def run_muscle(self):
        infile = os.path.join(self.wrkdir, self.id + ".pep")
        alnfile = infile + ".fasta"
        cmd = ["muscle", "-seqtype", "protein", "-quiet", 
                "-in", infile, "-out", alnfile]
        out = sp.run(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        logging.info(out.stderr.decode())
        #self.prot_aln = AlignIO.read(alnfile, "fasta")

    def run_trimal(self):
        alnfile = os.path.join(self.wrkdir, self.id + ".pep.fasta")
        cdsfile = os.path.join(self.wrkdir, self.id + ".cds")
        outfile = os.path.join(self.wrkdir, self.id + ".cds.fasta.trimal")
        cmd = ["trimal", "-in", alnfile, "-backtrans", cdsfile,
                "-out", outfile, "-fasta", "-automated1"]
        #print(' '.join(cmd))
        out = sp.run(cmd, stdout = sp.PIPE, stderr = sp.PIPE)
        logging.info(out.stderr.decode())
    
    def build_phylo(self):
        if self.phylotool == "raxml":
            self.run_raxml()
            #print(self.phylo)
            #Phylo.draw_ascii(self.phylo)
        elif self.phylotool == "phyml":
            pass
    
    def run_raxml(self):
        wd = os.path.join(os.path.abspath('.'), self.wrkdir)
        aln = os.path.join(wd, self.id + ".cds.fasta.trimal")
        #cmd = ["raxmlHPC-PTHREADS", "-T", "4", "-f", "a",
        cmd = ["raxmlHPC-PTHREADS", "-T", "1", "-f", "a",
                "-x", "601376", "-p", "601376", "-#", "100", "-w", wd,
        #        "-x", "601376", "-p", "601376", "-#", "10", "-w", wd,
                "-m", "GTRGAMMA", "-s", aln, "-n", self.id]
        logging.info(' '.join(cmd))
        out = sp.run(cmd, stderr = sp.PIPE, stdout = sp.PIPE)
        logging.info(out.stderr.decode())
        treefile = os.path.join(wd, 'RAxML_bipartitions.' + self.id)
        self.phylo = Phylo.read(treefile, 'newick')

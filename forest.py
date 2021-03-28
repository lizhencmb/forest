#!/usr/bin/env python
import click
import logging

@click.group()
def cli():
    pass

@cli.command()
@click.argument('genefamily', type=click.Path(exists=True))
@click.option('--protein', '-p', help='Protein file in fasta.')
@click.option('--cds', '-c', help='CDS file in fasta.')
@click.option('--wkdir', '-w', default='temp', help='Working directory.')
def trees(**kwargs):
    """Build gene trees"""
    _trees(**kwargs)

def _trees(genefamily, protein, cds, wkdir):
    from tools.genefamily import read_gene_families
    gfs = read_gene_families(genefamily, protein, cds, wkdir)
    for gf in gfs:
        gf.obtain_msa()
        gf.build_phylo()

@cli.command()
@click.argument('genefamily', type=click.Path(exists=True))
@click.option('-max', '-x', default=200, help='Maximum number of genes in a gene family.')
@click.option('-min', '-m', default=4, help='Minimum number of genes in a gene family.')
@click.option('-num', '-n', default=2, help='Number of split sets.')

def splitgf(**kwargs):
    """Split a large set of gene families"""
    _splitgf(**kwargs)

def _splitgf(genefamily, max, min, num):
    from tools.genefamily import split_gene_families
    split_gene_families(genefamily, max, min, num)

if __name__ == '__main__':
    cli()

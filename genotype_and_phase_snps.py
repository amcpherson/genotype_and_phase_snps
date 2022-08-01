import os
import click
import subprocess
from multiprocessing import Pool
from functools import partial


grch38_1kg_chromosomes = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
    'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']

grch38_1kg_vcf_filename_template = 'CCDG_14151_B01_GRM_WGS_2020-08-05_{chromosome}.filtered.shapeit2-duohmm-phased.vcf.gz'
grch38_1kg_bcf_filename_template = 'CCDG_14151_B01_GRM_WGS_2020-08-05_{chromosome}.filtered.shapeit2-duohmm-phased.bcf'
grch38_1kg_X_vcf_filename_template = 'CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz'
grch38_1kg_X_bcf_filename_template = 'CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.bcf'


def genotype_and_phase_chromosome(input_bam, input_reference, reference_dir, temps_prefix, chromosome):

    if chromosome == 'chrX':
        reference_vcf = os.path.join(reference_dir, grch38_1kg_X_vcf_filename_template)
        reference_bcf = os.path.join(reference_dir, grch38_1kg_X_bcf_filename_template)

    else:
        reference_vcf = os.path.join(reference_dir, grch38_1kg_vcf_filename_template.format(chromosome=chromosome))
        reference_bcf = os.path.join(reference_dir, grch38_1kg_bcf_filename_template.format(chromosome=chromosome))

    chromosome_mpileup_vcf = temps_prefix + '.' + chromosome + '.mpileup.vcf.gz'
    chromosome_calls_bcf = temps_prefix + '.' + chromosome + '.calls.bcf'

    subprocess.check_call([
        'bcftools', 'mpileup', '-Oz',
        '-f', input_reference,
        '--regions-file', reference_vcf,
        input_bam,
        '-o', chromosome_mpileup_vcf,
    ])

    subprocess.check_call([
        'bcftools', 'call', '-Ob', '-c',
        chromosome_mpileup_vcf,
        '-o', chromosome_calls_bcf,
    ])

    subprocess.check_call([
        'bcftools', 'index', chromosome_calls_bcf,
    ])

    genetic_map = os.path.join(reference_dir, f'{chromosome}.b38.gmap.gz')

    chromosome_phased_bcf = temps_prefix + '.' + chromosome + '.bcf'

    subprocess.check_call([
        'singularity', 'run', '--bind', '/juno', 'shapeit4_latest.sif',
        'shapeit4',
        '--input', chromosome_calls_bcf,
        '--map', genetic_map,
        '--region', chromosome,
        '--reference', reference_bcf,
        '--output', chromosome_phased_bcf,
        '--seed', '2'])

    return chromosome_phased_bcf


@click.command()
@click.argument('input_bam')
@click.argument('input_reference')
@click.argument('output_bcf')
@click.argument('reference_dir')
@click.argument('temps_prefix')
def main(input_bam, input_reference, output_bcf, reference_dir, temps_prefix):

    pool = Pool(len(grch38_1kg_chromosomes))
    phased_bcfs = pool.map(partial(genotype_and_phase_chromosome, input_bam, input_reference, reference_dir, temps_prefix), grch38_1kg_chromosomes)
    pool.close()
    pool.join()

    subprocess.check_call([
        'bcftools', 'concat', '-o', output_bcf] + phased_bcfs)
    

if __name__ == '__main__':
    main()


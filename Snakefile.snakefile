
from os.path import basename, dirname, join
from os import makedirs
import numpy as np

fpath_cna_input = config['cna_input']  # 23=X, 24=Y
fpath_cna_output = config['cna_output']
fpath_genomic_loc_by_varbin = config['genomic_loc_by_varbin']  #[0-based, 1-based]
fpath_genomic_loc_by_geneid = config['genomic_loc_by_geneid']  #[1-based, 1-based]

def _empty_str_2_nan(x):
    if x is None or x == '':
        return(np.nan)
    return(np.float(x))


rule all:
    input:
        fpath_cna_output,

rule _geneid_loc_1to0_based:
    input:
        txt = fpath_genomic_loc_by_geneid,
    output:
        unsorted_bed = temp(
            join('lib', basename(fpath_genomic_loc_by_geneid) + '.unsorted')),
        bed = join('lib', basename(fpath_genomic_loc_by_geneid) + '.bed'),
    run:
        fh_in = open(input.txt, 'r')
        fh_out = open(output.unsorted_bed, 'w+')
        next(fh_in)  # remove header
        for fin in fh_in:
            chrn, s1, e1, gid, gname, strd = fin.strip().split('\t')
            s0, e0 = int(s1) - 1, int(e1)
            if strd == '1':
                strd = '+'
            elif strd == '-1':
                strd = '-'
            else:
                strd = strd

            fh_out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                chrn, str(s0), str(e0),
                gid, gname, strd))
        fh_in.close()
        fh_out.close()
        shell('sort-bed {output.unsorted_bed} > {output.bed} ')

rule _varbin_loc_mixto0_based:
    input:
        txt = fpath_genomic_loc_by_varbin,
    output:
        unsorted_bed = temp(
            join('lib', basename(fpath_genomic_loc_by_varbin) + '.unsorted')),
        bed = join('lib', basename(fpath_genomic_loc_by_varbin) + '.bed3'),
        biosorted_bed = temp(
            join('lib', basename(fpath_genomic_loc_by_varbin) + '.sorted.bed3')),
    run:
        fh_in = open(input.txt, 'r')
        fh_out = open(output.unsorted_bed, 'w+')
        next(fh_in)  # remove header
        for fin in fh_in:
            chrn, s0, e1, _ = fin.strip().split('\t', 3)
            chrn = chrn.replace('chr', '')
            chrn = chrn.replace('Chr', '')

            e0 = int(e1) + 1
            fh_out.write('{}\t{}\t{}\n'.format(chrn, s0, e0))
        fh_in.close()
        fh_out.close()
        shell('sort-bed {output.unsorted_bed} > {output.bed} ')
        shell('sort -V -k1,1 -k2,2n -k3,3n {output.unsorted_bed} > {output.biosorted_bed}')

# Merge varbin with the CNA matrix
# shared key = varbin's chr + end = CNA's chr + chrompos
rule _assign_varbin_with_cna:
    input:
        varbin = rules._varbin_loc_mixto0_based.output.biosorted_bed,
        cna = fpath_cna_input,
    output:
        tsv = join(
            dirname(fpath_cna_output),
            basename(fpath_cna_input) + '.full_varbin'),
    run:
        fh_out = open(output.tsv, 'w+')

        dict_cna = {}
        dict_abs_pos = {}
        n_screens = -1
        with open(input.cna, 'r') as fr:
            fr_header = next(fr).strip('\n').split('\t')
            n_screens = len(fr_header) - 3  # exclude first three cols
            print('>> {:,} screens/cells/samples are available.'.format(n_screens))
            assert n_screens != -1, 'Expect >= 1 screens but here {}.'.format(n_screens)
            names_screens = '\t'.join(fr_header[3:])
            fh_out.write('chrom\tstart\tend\tabspos\t{}\n'.format(names_screens))

            for line in fr:
                chrn, chr_pos, g_pos, vals = line.strip('\n').split('\t', 3)
                if chrn == '23':
                    chrn = 'X'
                if chrn == '24':
                    chrn = 'Y'
                k_query = '{}\t{}'.format(chrn, chr_pos)
                dict_cna.setdefault(k_query, vals)
                dict_abs_pos.setdefault(k_query, g_pos)

        default_varbin_valinfo = '\t'.join(['' for i in range(n_screens)])
        with open(input.varbin, 'r') as fl:
            for line in fl:
                chrn, s, e = line.strip('\n').split('\t')
                k_query = '{}\t{}'.format(chrn, e)
                cna_query = dict_cna.get(k_query, default_varbin_valinfo)
                abspos_query = dict_abs_pos.get(k_query, 0)
                fh_out.write('{}\t{}\t{}\t{}\t{}\n'.format(
                    chrn, s, e, abspos_query, cna_query))
        fh_out.close()

rule _associate_geneid_to_varbin:
    input:
        varbin = rules._varbin_loc_mixto0_based.output.bed,
        geneid = rules._geneid_loc_1to0_based.output.bed,
    output:
        unsorted_txt = temp(join('lib', '{}-{}-{}.txt.unsorted'.format('associate', 'geneid', 'varbin'))),
        txt = join('lib', '{}-{}-{}.txt'.format('associate', 'geneid', 'varbin')),
    run:
        cmd = 'bedmap --delim \"\\t\" --echo --echo-map '
        cmd += '{input.geneid} {input.varbin} '
        # a gene is excluded if it overlaps with >=1 varbins
        # cmd += '| grep -v {} - '.format('";"')

        # a gene is excluded if it overlaps non varbins
        cmd += '| awk {} - '.format(' \'$7!=""\' ')
        cmd += '> {output.unsorted_txt} '
        shell(cmd)
        shell('sort -V -k1,1 -k2,2n -k3,3n {output.unsorted_txt} > {output.txt} ')


rule _cna_loc2geneid:
    input:
        gid_varbin = rules._associate_geneid_to_varbin.output.txt,
        varbin_cna = rules._assign_varbin_with_cna.output.tsv,
    output:
        txt = fpath_cna_output,
    run:
        import numpy as np
        fh_out = open(output.txt, 'w+')

        dict_varbin_cna = {}
        default_vals = None
        n_screens = -1
        with open(input.varbin_cna, 'r') as fh_in:
            header = next(fh_in).strip().split('\t')
            names_screens = header[4:]
            n_screens = len(names_screens)
            default_vals = np.array([np.nan for i in range(n_screens)])
            for line in fh_in:
                chrn, s, e, abspos, vals_str= line.strip('\n').split('\t', 4)
                vals_str = vals_str.strip('\n').split('\t')
                assert len(vals_str) == n_screens, 'parsing error.'

                vals = np.array(list(
                    map(lambda c: _empty_str_2_nan(c), vals_str)))
                k_query = '{}\t{}\t{}'.format(chrn, s, e)
                dict_varbin_cna.setdefault(k_query, vals)

        fh_out.write(
            'chrom\tstart\tend\tGENE_ID\tGENE_NAME\tstrand\t{}\n'.format(
                '\t'.join(names_screens)))

        with open(input.gid_varbin, 'r') as fh_in:
            for line in fh_in:
                g_chrn, g_s, g_e, gid, gname, gstrd, varbins = line.strip('\n').split('\t', 6)
                varbins_l = varbins.strip('\n').split(';')
                n_varbins = len(varbins_l)
                v_varbins = np.vstack([default_vals for i in range(n_varbins)])
                for i in range(n_varbins):
                    k_query = varbins_l[i]
                    v_query = dict_varbin_cna.get(k_query, default_vals)
                    v_varbins[i, :] = v_query

                v = np.array([np.nanmedian(v_varbins[:, j]) for j in range(n_screens)])
                v_str = '\t'.join(map(str, v)).replace('nan', '')
                fh_out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    g_chrn, g_s, g_e, gid, gname, gstrd, v_str))
        fh_out.close()

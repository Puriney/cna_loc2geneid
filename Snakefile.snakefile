
from os.path import basename, dirname, join
from os import makedirs

fpath_cna_input = config['cna_input']  # 23=X, 24=Y
fpath_cna_output = config['cna_output']
fpath_genomic_loc_by_varbin = config['genomic_loc_by_varbin']  #[0-based, 1-based]
fpath_genomic_loc_by_geneid = config['genomic_loc_by_geneid']  #[1-based, 1-based]

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
        bed = join('lib', basename(fpath_genomic_loc_by_varbin) + '.bed3'),
    run:
        fh_in = open(input.txt, 'r')
        fh_out = open(output.bed, 'w+')
        next(fh_in)  # remove header
        for fin in fh_in:
            chrn, s0, e1, _ = fin.strip().split('\t', 3)
            chrn = chrn.replace('chr', '')
            chrn = chrn.replace('Chr', '')

            e0 = int(e1) + 1
            fh_out.write('{}\t{}\t{}\n'.format(chrn, s0, e0))
        fh_in.close()
        fh_out.close()

# Merge varbin with the CNA matrix
# shared key = varbin's chr + end = CNA's chr + chrompos
rule _assign_varbin_with_cna:
    input:
        varbin = rules._varbin_loc_mixto0_based.output.bed,
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
        txt = join('lib', '{}-{}-{}.txt'.format('associate', 'geneid', 'varbin')),
    run:
        cmd = 'bedmap --delim \"\\t\" --echo --echo-map '
        cmd += '{input.geneid} {input.varbin} '
        cmd += '| grep -v {} - '.format('";"')
        cmd += '| awk {} - '.format(' \'$7!=""\' ')
        cmd += '> {output.txt} '
        shell(cmd)

rule _cna_loc2geneid:
    input:
        gid_varbin = rules._associate_geneid_to_varbin.output.txt,
        varbin_cna = rules._assign_varbin_with_cna.output.tsv,
    output:
        txt = fpath_cna_output,
    run:
        fh_out = open(output.txt, 'w+')

        dict_varbin_cna = {}
        default_vals = ''
        with open(input.varbin_cna, 'r') as fh_in:
            header = next(fh_in).strip().split('\t')
            names_screens = header[4:]
            default_vals = '\t'.join(['' for i in range(len(names_screens))])
            for line in fh_in:
                chrn, s, e, abspos, vals= line.strip('\n').split('\t', 4)
                k_query = '{}\t{}\t{}'.format(chrn, s, e)
                dict_varbin_cna.setdefault(k_query, vals)

        fh_out.write(
            'chrom\tstart\tend\tGENE_ID\tGENE_NAME\tstrand\t{}\n'.format(
                '\t'.join(names_screens)))

        with open(input.gid_varbin, 'r') as fh_in:
            for line in fh_in:
                g_chrn, g_s, g_e, gid, gname, gstrd, varbin = line.strip().split('\t', 6)
                k_query = varbin
                v_query = dict_varbin_cna.get(k_query, default_vals)
                fh_out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    g_chrn, g_s, g_e, gid, gname, gstrd, v_query))
        fh_out.close()


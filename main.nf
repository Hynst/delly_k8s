process delly_somatic {

    publishDir "${params.pubdir}/results/SV_delly", mode: 'copy'
    container "dellytools/delly:latest"

    input:
        tuple val(sample), val(t_bam), val(t_bai), val(n_bam), val(n_bai)
    
    output:
        path '.'

    // sample here is a string -> tumorID_normalID
    shell:
    """
    delly call -o ${sample}_SV_delly.bcf -g ${params.ref} ${t_bam} ${n_bam}
    """

}

workflow {
    // create channel
    input_ch = Channel.empty()
    tsv = file(params.input)
    input_ch = extractFastq(tsv)

    delly_somatic(input_ch)

}


def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
}

def extractFastq(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->        
            def sample      = row[0]
            def t_bam       = returnFile(row[1])
            def t_bai       = returnFile(row[2])
            def n_bam       = returnFile(row[3])
            def n_bai       = returnFile(row[4])

            [sample, t_bam, t_bai, n_bam, n_bai]
        }
}

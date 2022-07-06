task participant_vcfs {

    File vcf_file
    String participant_id

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        date +"[%b %d %H:%M:%S] Generating participant VCF (SNPs only)"
        # select SNPs, filter out missing sites
        bcftools view --no-update -s ${participant_id} -v snps ${vcf_file} | bcftools view --no-update -e 'GT=".|."' -Oz -o ${participant_id}.snps.vcf.gz
        tabix ${participant_id}.snps.vcf.gz

        date +"[%b %d %H:%M:%S] Subsetting biallelic het sites for ASE"
        bcftools view --no-update -i 'GT="het"' ${participant_id}.snps.vcf.gz | bcftools norm -m+ | bcftools view -m2 -M2 -Oz -o ${participant_id}.snps.het.vcf.gz
        tabix ${participant_id}.snps.het.vcf.gz

        date +"[%b %d %H:%M:%S] Done"
    }

    output {
        File snps_vcf = "${participant_id}.snps.vcf.gz"
        File snps_vcf_index = "${participant_id}.snps.vcf.gz.tbi"
        File snps_het_vcf = "${participant_id}.snps.het.vcf.gz"
        File snps_het_vcf_index = "${participant_id}.snps.het.vcf.gz.tbi"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V10"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow participant_vcfs_workflow {
    call participant_vcfs
}

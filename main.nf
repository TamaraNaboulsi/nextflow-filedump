#!/usr/bin/env nextflow

nextflow.enable.dsl=2

workflow {
  DbFactory()

  if (params.run_datachecks) {
    RunDataChecks(DbFactory.out.splitText())
    if (params.dump_mysql) {
      dump_mysql(RunDataChecks.out.splitText())
      SpeciesFactory(RunDataChecks.out)
    } else {
      SpeciesFactory(RunDataChecks.out)
    }
  } else if (params.dump_mysql && !params.run_datachecks) {
    dump_mysql(DbFactory.out.splitText())
    SpeciesFactory(DbFactory.out.splitText())
  } else {
    SpeciesFactory(DbFactory.out.splitText())
  }

  genome(SpeciesFactory.out.splitText())
  geneset(SpeciesFactory.out.splitText())
  rna_seq(SpeciesFactory.out.splitText())

  sync_data = channel.empty()
  if (params.dump_mysql) {
    sync_data = sync_data.concat(dump_mysql.out, genome.out, geneset.out, rna_seq.out)
  } else {
    sync_data = sync_data.concat(genome.out, geneset.out, rna_seq.out)
  }

  sync(sync_data)

  Metadata_JSON(sync.out)
  if (params.ftp_root) {
    Sync_Metadata(Metadata_JSON.out.splitText())
  }
}

workflow dump_mysql {
  take: dump_mysql_in
  main:
    MySQL_TXT(dump_mysql_in.splitText())
    MySQL_Compress(MySQL_TXT.out[0].splitText())
    Checksum(MySQL_TXT.out[1].splitText(), MySQL_Compress.out.collectFile())
    Verify(Checksum.out.splitText())
  emit:
    dump_mysql_out = Verify.out.collectFile()
}

workflow genome {
  take: genome_in
  main:
    GenomeDirectoryPaths(genome_in.splitText())

    genome_data = channel.empty()
    if (params.run_assembly_chain) {
      Assembly_Chain(GenomeDirectoryPaths.out.splitText())
      genome_data = genome_data.concat(Assembly_Chain.out)
    }
    if (params.run_chromosome_tsv) {
      Chromosome_TSV(GenomeDirectoryPaths.out.splitText())
      genome_data = genome_data.concat(Chromosome_TSV.out)
    }
    if (params.run_genome_fasta) {
      Genome_FASTA(GenomeDirectoryPaths.out.splitText())
      genome_data = genome_data.concat(Genome_FASTA.out[0])
    }

    Genome_Compress(genome_data.splitText())

    if (params.run_genome_fasta) {
      Symlink_Genome_FASTA(Genome_FASTA.out[0].splitText(), Genome_FASTA.out[1], Genome_FASTA.out[2], Genome_Compress.out.collectFile())
      Checksum(GenomeDirectoryPaths.out.splitText(), Symlink_Genome_FASTA.out.collectFile())
    } else {
      Checksum(GenomeDirectoryPaths.out.splitText(), Genome_Compress.out.collectFile())
    }
    Verify(Checksum.out.splitText())
  emit:
    genome_out = Verify.out.collectFile()
}

workflow geneset {
  take: geneset_in
  main:
    GenesetDirectoryPaths(geneset_in.splitText())

    geneset_data = channel.empty()
    if (params.run_geneset_embl) {
      Geneset_EMBL(GenesetDirectoryPaths.out.splitText())
      geneset_data = geneset_data.concat(Geneset_EMBL.out)
    }
    if (params.run_geneset_fasta) {
      Geneset_FASTA(GenesetDirectoryPaths.out.splitText())
      geneset_data = geneset_data.concat(Geneset_FASTA.out)
    }
    if (params.run_geneset_gff3) {
      Geneset_GFF3(GenesetDirectoryPaths.out.splitText())
      geneset_data = geneset_data.concat(Geneset_GFF3.out)
    }
    if (params.run_geneset_gtf) {
      Geneset_GTF(GenesetDirectoryPaths.out.splitText())
      geneset_data = geneset_data.concat(Geneset_GTF.out)
    }
    if (params.run_xref_tsv) {
      Xref_TSV(GenesetDirectoryPaths.out.splitText())
      geneset_data = geneset_data.concat(Xref_TSV.out)
    }

    Geneset_Compress(geneset_data.splitText())

    if (params.run_xref_tsv) {
      Symlink_Xref_TSV(Xref_TSV.out[0].splitText(), Xref_TSV.out[1], Xref_TSV.out[2], Geneset_Compress.out.collectFile())
      Checksum(GenesetDirectoryPaths.out.splitText(), Symlink_Xref_TSV.out.collectFile())
    } else {
      Checksum(GenesetDirectoryPaths.out.splitText(), Geneset_Compress.out.collectFile())
    }
    Verify(Checksum.out.splitText())
  emit:
    geneset_out = Verify.out.collectFile()
}

workflow rna_seq {
  take: rna_seq_in
  main:
    RNASeqDirectoryPaths(rna_seq_in.splitText())

    if (params.run_rnaseq_exists) {
      RNASeq_Exists(RNASeqDirectoryPaths.out.splitText())
      if (RNASeq_Exists.out[1] == 'exists') {
        Symlink_RNASeq(RNASeq_Exists.out[0].splitText())
        Verify_Unzipped(Symlink_RNASeq.out.collectFile())
        RNASeq_Missing(RNASeq_Exists.out[2].splitText())
        rna_seq_out = Verify_Unzipped.out.collectFile()
      } else {
        rna_seq_out = channel.empty()
      }
    }
  emit:
    rna_seq_out
}

workflow sync {
  take: sync_in
  main:
    Sync(sync_in.splitText())
    README(Sync.out.splitText())
  emit:
    sync_out = README.out.collectFile()
}

process DbFactory {
  output:
  path 'dataflow_2.json'

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::DbFactory' -reg_conf=${params.registry} -species=${params.species} -antispecies=${params.antispecies} -division=${params.division} -run_all=${params.run_all} -dbname=${params.dbname} -meta_filters=${params.meta_filters}
  """
}

process RunDataChecks {
  input:
  each x

  output:
  val x

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::RunDataChecks' -dataflow='$x' -registry_file=${params.registry} -history_file=${params.history_file} -output_dir=${params.output_dir} -config_file=${params.config_file} -datacheck_names=${params.datacheck_names} -datacheck_groups=${params.datacheck_groups} -datacheck_types=${params.datacheck_types} -failures_fatal=1
  """
}

process MySQL_TXT {
  input:
  each x

  output:
  path 'dataflow_2.json'
  path 'dataflow_3.json'

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::MySQL_TXT' -dataflow='$x' -reg_conf=${params.registry} -dump_dir=${params.dump_dir} -overwrite=${params.overwrite} -output_dir=${params.output_dir}
  """
}

process SpeciesFactory {
  input:
  each x

  output:
  path 'dataflow_2.json'

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::DbAwareSpeciesFactory' -dataflow='$x' -reg_conf=${params.registry}
  """
}

process MySQL_Compress {
  input:
  each x

  output:
  val 'done'

  shell:
  output_filename = (x =~ /"output_filename":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  script:
  """
  gzip -n -f $output_filename
  """
}

process Checksum {
  input:
  each x
  val y

  output:
  val x

  shell:
  output_dir = (x =~ /"output_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  """
  cd $output_dir; find -L . -type f -not -name "md5sum.txt" | sed 's#^./##' | xargs md5sum > md5sum.txt
  """
}

process Verify {
  input:
  each x

  output:
  val x

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Verify' -dataflow='$x'
  """
}

process Sync {
  input:
  each x

  output:
  val x

  shell:
  output_dir = (x =~ /"output_dir":"([A-Za-z0-9_.-\/]+)"/)[0][1]
  if (x =~ /"ftp_dir"/)
    ftp_dir = (x =~ /"ftp_dir":"*([A-Za-z0-9_.\-\/]+)"*/)[0][1]

  """
  if [ -z "$ftp_dir" ]; then
    mkdir -p $ftp_dir; rsync -aLW $output_dir/ $ftp_dir
  fi
  """
}

process README {
  input:
  each x

  output:
  val 'done'

  shell:
  data_category = (x =~ /"data_category":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  script:
  if (data_category == 'geneset' || data_category == 'genome')
    """
    perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::README' -dataflow='$x'
    """
}

// Genome

process GenomeDirectoryPaths {
  input:
  each x

  output:
  path 'dataflow_3.json'

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::DirectoryPaths' -dataflow='$x' -reg_conf=${params.registry} -data_category='genome' -analysis_types=${params.genome_types} -dump_dir=${params.dump_dir}
  """
}

process Assembly_Chain {
  input:
  each x

  output:
  path 'dataflow_2.json'

  shell:
  species_name = (x =~ /"species_name":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Assembly_Chain' -dataflow='$x' -reg_conf=${params.registry} -ucsc=${params.chain_ucsc} -overwrite=${params.overwrite} -species=$species_name
  """
}

process Chromosome_TSV {
  input:
  each x

  output:
  path 'dataflow_2.json'

  shell:
  species_name = (x =~ /"species_name":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Chromosome_TSV' -dataflow='$x' -reg_conf=${params.registry} -overwrite=${params.overwrite} -species=$species_name
  """
}

process Genome_FASTA {
  memory { 1.GB * task.attempt }
  errorStrategy 'retry'

  input:
  each x

  output:
  path 'dataflow_2.json'
  val output_dir
  val data_category

  shell:
  output_dir = (x =~ /"output_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  species_name = (x =~ /"species_name":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  data_category = (x =~ /"data_category":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Genome_FASTA' -dataflow='$x' -reg_conf=${params.registry} -blast_index=${params.blast_index} -blastdb_exe=${params.blastdb_exe} -per_chromosome=${params.per_chromosome} -overwrite=${params.overwrite} -species=$species_name
  """
}

process Genome_Compress {
  memory { 1.GB * task.attempt }
  errorStrategy 'retry'

  input:
  each x

  output:
  val x

  shell:
  output_filename = (x =~ /"output_filename":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  """
  if [ -f "$output_filename" ]; then
    gzip -n -f $output_filename
  fi
  """
}

process Symlink_Genome_FASTA {
  input:
  each x
  val output_dir
  val data_category
  val y

  output:
  val 'done'

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Symlink' -dataflow='$x' -data_category=$data_category -dump_dir=${params.dump_dir} -output_dir=$output_dir
  """
}

// Geneset

process GenesetDirectoryPaths {
  input:
  each x

  output:
  path 'dataflow_3.json'

  script:
  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::DirectoryPaths' -dataflow='$x' -reg_conf=${params.registry} -data_category='geneset' -analysis_types=${params.geneset_types} -dump_dir=${params.dump_dir}
  """
}

process Geneset_EMBL {
  memory { 1.GB * task.attempt }
  errorStrategy 'retry'

  input:
  each x

  output:
  path 'dataflow_2.json'

  shell:
  species_name = (x =~ /"species_name":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Geneset_EMBL' -dataflow='$x' -reg_conf=${params.registry} -per_chromosome=${params.embl_per_chromosome} -overwrite=${params.overwrite} -species=$species_name
  """
}

process Geneset_FASTA {
  memory { 1.GB * task.attempt }
  errorStrategy 'retry'

  input:
  each x

  output:
  path 'dataflow_2.json'

  shell:
  species_name = (x =~ /"species_name":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Geneset_FASTA' -dataflow='$x' -reg_conf=${params.registry} -blast_index=${params.blast_index} -blastdb_exe=${params.blastdb_exe} -overwrite=${params.overwrite} -species=$species_name
  """
}

process Geneset_GFF3 {
  memory { 1.GB * task.attempt }
  errorStrategy 'retry'

  input:
  each x

  output:
  path 'dataflow_2.json'

  shell:
  species_name = (x =~ /"species_name":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Geneset_GFF3' -dataflow='$x' -reg_conf=${params.registry} -per_chromosome=${params.gff3_per_chromosome} -gt_gff3_exe=${params.gt_gff3_exe} -gt_gff3validator_exe=${params.gt_gff3validator_exe} -overwrite=${params.overwrite} -species=$species_name
  """
}

process Geneset_GTF {
  memory { 1.GB * task.attempt }
  errorStrategy 'retry'

  input:
  each x

  output:
  path 'dataflow_2.json'

  shell:
  species_name = (x =~ /"species_name":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Geneset_GTF' -dataflow='$x' -reg_conf=${params.registry} -per_chromosome=${params.gtf_per_chromosome} -gtf_to_genepred_exe=${params.gtf_to_genepred_exe} -genepred_check_exe=${params.genepred_check_exe} -overwrite=${params.overwrite} -species=$species_name
  """
}

process Xref_TSV {
  memory { 1.GB * task.attempt }
  errorStrategy 'retry'

  input:
  each x

  output:
  path 'dataflow_2.json'
  val output_dir
  val data_category

  shell:
  output_dir = (x =~ /"output_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  species_name = (x =~ /"species_name":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  data_category = (x =~ /"data_category":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Xref_TSV' -dataflow='$x' -reg_conf=${params.registry} -external_dbs=${params.xref_external_dbs} -overwrite=${params.overwrite} -species=$species_name
  """
}

process Geneset_Compress {
  memory { 1.GB * task.attempt }
  errorStrategy 'retry'

  input:
  each x

  output:
  val x

  shell:
  output_filename = (x =~ /"output_filename":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  """
  gzip -n -f $output_filename
  """
}

process Symlink_Xref_TSV {
  input:
  each x
  val output_dir
  val data_category
  val y

  output:
  val 'done'

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Symlink' -dataflow='$x' -data_category=$data_category -dump_dir=${params.dump_dir} -output_dir=$output_dir
  """
}

// RNASeq

process RNASeqDirectoryPaths {
  input:
  each x

  output:
  path 'dataflow_3.json'

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::DirectoryPaths' -dataflow='$x' -reg_conf=${params.registry} -data_category='rnaseq' -analysis_types=${params.rnaseq_types} -dump_dir=${params.dump_dir}
  """
}

process RNASeq_Exists {
  input:
  each x

  output:
  val x
  val exists
  path 'dataflow_3.json'

  shell:
  species_name = (x =~ /"species_name":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  exists = 'no'

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::RNASeq_Exists' -dataflow='$x' -reg_conf=${params.registry} -overwrite=${params.overwrite} -species=$species_name

  if [ -f "dataflow_2.json" ]; then
    $exists='yes'
  fi
  if [ ! -f "dataflow_3.json" ]; then
    echo '' > dataflow_3.json
  fi
  """
}

process Symlink_RNASeq {
  input:
  each x

  output:
  val x

  shell:
  species_name = (x =~ /"species_name":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Symlink_RNASeq' -dataflow='$x' -reg_conf=${params.registry} -overwrite=${params.overwrite} -dump_dir=${params.dump_dir} -species=$species_name
  """
}

process Verify_Unzipped {
  input:
  each x

  output:
  val x

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Verify' -dataflow='$x' -check_unzipped=0
  """
}

process RNASeq_Missing {
  input:
  each x

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::RNASeq_Missing' -dataflow='$x'
  """
}

// Metadata

process Metadata_JSON {
  input:
  val x

  output:
  path 'dataflow_2.json'

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Metadata_JSON' -dump_dir=${params.dump_dir}
  """
}

process Sync_Metadata {
  input:
  each x

  shell:
  output_filename = (x =~ /"output_filename":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  """
  rsync -aLW $output_filename ${params.ftp_root}
  """
}

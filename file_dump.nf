#!/usr/bin/env nextflow

nextflow.enable.dsl=2

workflow {
  DbFactory()

  if (params.run_datachecks) {
    RunDataChecks(DbFactory.out.splitText())
    if (params.dump_mysql) {
      MySQL_TXT(RunDataChecks.out)
      SpeciesFactory(RunDataChecks.out)
    } else {
      SpeciesFactory(RunDataChecks.out)
    }
  } else if (params.dump_mysql && !params.run_datachecks) {
    MySQL_TXT(DbFactory.out.splitText())
    SpeciesFactory(DbFactory.out.splitText())
  } else {
    SpeciesFactory(DbFactory.out.splitText())
  }

  if (params.dump_mysql) {
    MySQL_Compress(MySQL_TXT.out[0].splitText())
    Checksum(MySQL_TXT.out[1].splitText(), MySQL_Compress.out.collectFile())
  }

  GenomeDirectoryPaths(SpeciesFactory.out.splitText())
  GenesetDirectoryPaths(SpeciesFactory.out.splitText())
  RNASeqDirectoryPaths(SpeciesFactory.out.splitText())

  checksum_data = channel.empty()

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
    genome_data = genome_data.concat(Genome_FASTA.out)
  }
  Genome_Compress(genome_data.splitText())
  if (Genome_Compress.out) {
    if (params.run_genome_fasta) {
      Symlink_Genome_FASTA(Genome_FASTA.out[0].splitText(), Genome_FASTA.out[1], Genome_FASTA.out[2])
      Symlink_Genome_FASTA.out.view()
      if (Symlink_Genome_FASTA.out) {
        checksum_data = checksum_data.concat(GenomeDirectoryPaths.out)
      }
    } else {
      checksum_data.concat(GenomeDirectoryPaths.out)
    }
  }
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

  shell:
  species = (x =~ /"species":"([A-Za-z0-9_.\-]+)"/)[0][1]
  dbname = (x =~ /"dbname":"([A-Za-z0-9_.\-]+)"/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::RunDataChecks' -registry_file=${params.registry} -history_file=${params.history_file} -output_dir=${params.output_dir} -config_file=${params.config_file} -datacheck_names=${params.datacheck_names} -datacheck_groups=${params.datacheck_groups} -datacheck_types=${params.datacheck_types} -failures_fatal=1 -species=$species -dbname=$dbname
  """
}

process MySQL_TXT {
  input:
  each x

  output:
  path 'dataflow_2.json'
  path 'dataflow_3.json'

  shell:
  species = (x =~ /"species":"([A-Za-z0-9_.\-]+)"/)[0][1]
  dbname = (x =~ /"dbname":"([A-Za-z0-9_.\-]+)"/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::MySQL_TXT' -reg_conf=${params.registry} -dump_dir=${params.dump_dir} -overwrite=${params.overwrite} -output_dir=${params.output_dir} -species=$species -dbname=$dbname
  """
}

process SpeciesFactory {
  input:
  each x

  output:
  path 'dataflow_2.json'

  shell:
  species_list = (x =~ /"species_list":\["([A-Za-z0-9_.\-\/]+)"\]/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::DbAwareSpeciesFactory' -reg_conf=${params.registry} -species_list=$species_list
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
  val x
  val y

  output:
  val x

  shell:
  output_dir = (x =~ /"output_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  """
  cd $output_dir; find -L . -type f -not -name "md5sum.txt" | sed 's#^./##' | xargs md5sum > md5sum.txt
  """
}

// Genome

process GenomeDirectoryPaths {
  input:
  each x

  output:
  path 'dataflow_3.json'

  shell:
  species = (x =~ /"species":"([A-Za-z0-9_.\-]+)"/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::DirectoryPaths' -reg_conf=${params.registry} -data_category='genome' -analysis_types=${params.genome_types} -species=$species -dump_dir=${params.dump_dir}
  """
}

process Assembly_Chain {
  input:
  each x

  output:
  path 'dataflow_2.json'

  shell:
  output_dir = (x =~ /"output_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  web_dir = (x =~ /"web_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  species_name = (x =~ /"species_name":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  assembly = (x =~ /"assembly":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  data_category = (x =~ /"data_category":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  timestamped_dir = (x =~ /"timestamped_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  ftp_dir = (x =~ /"ftp_dir":"*([A-Za-z0-9_.\-\/]+)"*/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Assembly_Chain' -reg_conf=${params.registry} -ucsc=${params.chain_ucsc} -overwrite=${params.overwrite} -output_dir=$output_dir -web_dir=$web_dir -species_name=$species_name -species=$species_name -assembly=$assembly -data_category=$data_category -timestamped_dir=$timestamped_dir -ftp_dir=$ftp_dir
  """
}

process Chromosome_TSV {
  input:
  each x

  output:
  path 'dataflow_2.json'

  shell:
  output_dir = (x =~ /"output_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  web_dir = (x =~ /"web_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  species_name = (x =~ /"species_name":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  assembly = (x =~ /"assembly":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  data_category = (x =~ /"data_category":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  timestamped_dir = (x =~ /"timestamped_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  ftp_dir = (x =~ /"ftp_dir":"*([A-Za-z0-9_.\-\/]+)"*/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Chromosome_TSV' -reg_conf=${params.registry} -overwrite=${params.overwrite} -output_dir=$output_dir -web_dir=$web_dir -species_name=$species_name -species=$species_name -assembly=$assembly -data_category=$data_category -timestamped_dir=$timestamped_dir -ftp_dir=$ftp_dir
  """
}

process Genome_FASTA {
  input:
  each x

  output:
  path 'dataflow_2.json'
  val output_dir
  val data_category

  shell:
  output_dir = (x =~ /"output_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  web_dir = (x =~ /"web_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  species_name = (x =~ /"species_name":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  assembly = (x =~ /"assembly":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  data_category = (x =~ /"data_category":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  timestamped_dir = (x =~ /"timestamped_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  ftp_dir = (x =~ /"ftp_dir":"*([A-Za-z0-9_.\-\/]+)"*/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Genome_FASTA' -reg_conf=${params.registry} -blast_index=${params.blast_index} -blastdb_exe=${params.blastdb_exe} -per_chromosome=${params.per_chromosome} -overwrite=1 -output_dir=$output_dir -web_dir=$web_dir -species_name=$species_name -species=$species_name -assembly=$assembly -data_category=$data_category -timestamped_dir=$timestamped_dir -ftp_dir=$ftp_dir
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
  val y
  val z

  output:
  val x

  shell:
  output_filename = (x =~ /"output_filename":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  """
  echo $x
  echo $y
  """
}

// Geneset

process GenesetDirectoryPaths {
  input:
  each x

  output:
  path 'dataflow_3.json'

  shell:
  species = (x =~ /"species":"([A-Za-z0-9_.\-]+)"/)[0][1]

  script:
  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::DirectoryPaths' -reg_conf=${params.registry} -data_category='geneset' -analysis_types=${params.geneset_types} -species=$species -dump_dir=${params.dump_dir}
  """
}

process Geneset_EMBL {
  input:
  each x

  output:
  path 'dataflow_2.json'

  shell:
  output_dir = (x =~ /"output_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  web_dir = (x =~ /"web_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  species_name = (x =~ /"species_name":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  assembly = (x =~ /"assembly":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  data_category = (x =~ /"data_category":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  timestamped_dir = (x =~ /"timestamped_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  ftp_dir = (x =~ /"ftp_dir":"*([A-Za-z0-9_.\-\/]+)"*/)[0][1]
  geneset = (x =~ /"geneset":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Geneset_EMBL' -reg_conf=${params.registry} -per_chromosome=${params.embl_per_chromosome} -overwrite=${params.overwrite} -output_dir=$output_dir -web_dir=$web_dir -species_name=$species_name -species=$species_name -assembly=$assembly -data_category=$data_category -timestamped_dir=$timestamped_dir -ftp_dir=$ftp_dir -geneset=$geneset
  """
}

process Geneset_FASTA {
  input:
  each x

  output:
  path 'dataflow_2.json'

  shell:
  output_dir = (x =~ /"output_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  web_dir = (x =~ /"web_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  species_name = (x =~ /"species_name":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  assembly = (x =~ /"assembly":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  data_category = (x =~ /"data_category":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  timestamped_dir = (x =~ /"timestamped_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  ftp_dir = (x =~ /"ftp_dir":"*([A-Za-z0-9_.\-\/]+)"*/)[0][1]
  geneset = (x =~ /"geneset":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Geneset_FASTA' -reg_conf=${params.registry} -blast_index=${params.blast_index} -blastdb_exe=${params.blastdb_exe} -overwrite=${params.overwrite} -output_dir=$output_dir -web_dir=$web_dir -species_name=$species_name -species=$species_name -assembly=$assembly -data_category=$data_category -timestamped_dir=$timestamped_dir -ftp_dir=$ftp_dir -geneset=$geneset
  """
}

process Geneset_GFF3 {
  input:
  each x

  output:
  path 'dataflow_2.json'

  shell:
  output_dir = (x =~ /"output_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  web_dir = (x =~ /"web_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  species_name = (x =~ /"species_name":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  assembly = (x =~ /"assembly":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  data_category = (x =~ /"data_category":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  timestamped_dir = (x =~ /"timestamped_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  ftp_dir = (x =~ /"ftp_dir":"*([A-Za-z0-9_.\-\/]+)"*/)[0][1]
  geneset = (x =~ /"geneset":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Geneset_GFF3' -reg_conf=${params.registry} -per_chromosome=${params.gff3_per_chromosome} -gt_gff3_exe=${params.gt_gff3_exe} -gt_gff3validator_exe=${params.gt_gff3validator_exe} -overwrite=${params.overwrite} -output_dir=$output_dir -web_dir=$web_dir -species_name=$species_name -species=$species_name -assembly=$assembly -data_category=$data_category -timestamped_dir=$timestamped_dir -ftp_dir=$ftp_dir -geneset=$geneset
  """
}

process Geneset_GTF {
  input:
  each x

  output:
  path 'dataflow_2.json'

  shell:
  output_dir = (x =~ /"output_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  web_dir = (x =~ /"web_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  species_name = (x =~ /"species_name":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  assembly = (x =~ /"assembly":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  data_category = (x =~ /"data_category":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  timestamped_dir = (x =~ /"timestamped_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  ftp_dir = (x =~ /"ftp_dir":"*([A-Za-z0-9_.\-\/]+)"*/)[0][1]
  geneset = (x =~ /"geneset":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Geneset_GTF' -reg_conf=${params.registry} -per_chromosome=${params.gtf_per_chromosome} -gtf_to_genepred_exe=${params.gtf_to_genepred_exe} -genepred_check_exe=${params.genepred_check_exe} -overwrite=${params.overwrite} -output_dir=$output_dir -web_dir=$web_dir -species_name=$species_name -species=$species_name -assembly=$assembly -data_category=$data_category -timestamped_dir=$timestamped_dir -ftp_dir=$ftp_dir -geneset=$geneset
  """
}

process Xref_TSV {
  input:
  each x

  output:
  path 'dataflow_2.json'

  shell:
  output_dir = (x =~ /"output_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  web_dir = (x =~ /"web_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  species_name = (x =~ /"species_name":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  assembly = (x =~ /"assembly":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  data_category = (x =~ /"data_category":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  timestamped_dir = (x =~ /"timestamped_dir":"([A-Za-z0-9_.\-\/]+)"/)[0][1]
  ftp_dir = (x =~ /"ftp_dir":"*([A-Za-z0-9_.\-\/]+)"*/)[0][1]
  geneset = (x =~ /"geneset":"([A-Za-z0-9_.\-\/]+)"/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::Xref_TSV' -reg_conf=${params.registry} -external_dbs=${params.xref_external_dbs} -overwrite=${params.overwrite} -output_dir=$output_dir -web_dir=$web_dir -species_name=$species_name -species=$species_name -assembly=$assembly -data_category=$data_category -timestamped_dir=$timestamped_dir -ftp_dir=$ftp_dir -geneset=$geneset
  """
}

// RNASeq

process RNASeqDirectoryPaths {
  input:
  each x

  output:
  val x

  shell:
  species = (x =~ /"species":"([A-Za-z0-9_.\-]+)"/)[0][1]

  """
  perl ${params.pipeline_dir}/run_process.pl -class='Nextflow::FileDump::DirectoryPaths' -reg_conf=${params.registry} -data_category='rnaseq' -analysis_types=${params.rnaseq_types} -species=$species -dump_dir=${params.dump_dir}
  """
}

process RNASeq_Exists {
  input:
  each x

  output:
  path 'dataflow_2.json'

  shell:


  """

  """
}

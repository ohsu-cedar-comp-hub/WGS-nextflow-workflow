#!/usr/bin/env nextflow

include {LearnReadOrientationModel} from '../tools/gatk/learn_read_orientation_model.nf'

//Define workflow for learn read orientation model

workflow { 

	LearnReadOrientationModel (file(params.f1r2_tar_gz), "test")

  }

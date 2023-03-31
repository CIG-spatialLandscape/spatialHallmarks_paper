OCT
spaceranger count --id=sample_id --transcriptome=../refdata-gex-GRCh38-2020-A --fastqs=directoy_to_fastq_files --sample=sample_name --image=tiff_image --slide=slide_name --slidefile=slide_file(.gpr) --area=A1

FFPE
spaceranger count --id=sample_id --transcriptome=../refdata-gex-GRCh38-2020-A --probe-set=../Visium_Human_Transcriptome_Probe_Set_v1$ --image=tiff_image --slide=slide_name --slide=slide_name --slidefile=slide_file(.gpr) --area=B1 --reorient-images
#--reorient-images only if image is not properly oriented

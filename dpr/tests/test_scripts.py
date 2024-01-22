import subprocess
import pytest

from pathlib import Path

cwd = Path(__file__).parents[2] / Path("datasets")
commands_dpr = [
    "dpr data_SENSE3_MB3_noisemap.nii.gz sigma_maxlk_nmaps.nii.gz N_maxlk_nmaps.nii.gz mask_maxlk_nmaps.nii.gz -m maxlk --noise_maps",
    "dpr data_SENSE3_MB3_noisemap.nii.gz sigma_nmaps.nii.gz N_nmaps.nii.gz mask_nmaps.nii.gz --noise_maps",
    "dpr data_SENSE3_MB3_noisemap.nii.gz sigma_nmaps.nii.gz N_nmaps.nii.gz mask_nmaps.nii.gz --noise_maps -f --subsample",
    "dpr data_SENSE3_MB3_noisemap.nii.gz sigma_nmaps.nii.gz N_nmaps.nii.gz mask_nmaps.nii.gz --noise_maps -f --fast_median -m maxlk",
    "dpr dwi_1_8.nii.gz sigma.nii.gz N.nii.gz mask.nii.gz -v",
    "dpr dwi_1_8.nii.gz sigma.nii.gz N.nii.gz mask.nii.gz -m maxlk -f --ncores 4",
    "dpr dwi_1_8.nii.gz sigma_maxlk.nii.gz N_maxlk.nii.gz mask_maxlk.nii.gz -m maxlk --size 3 -f -v --axis 0",
]

commands_dpr_graph = [
    "dpr data_SENSE3_MB3_noisemap.nii.gz sigma_maxlk_nmaps.nii.gz N_maxlk_nmaps.nii.gz mask_maxlk_nmaps.nii.gz -m maxlk --noise_maps",
    "dpr data_SENSE3_MB3_noisemap.nii.gz sigma_nmaps.nii.gz N_nmaps.nii.gz mask_nmaps.nii.gz --noise_maps",
    "dpr data_SENSE3_MB3_noisemap.nii.gz sigma_nmaps.nii.gz N_nmaps.nii.gz mask_nmaps.nii.gz --noise_maps -f --subsample",
    "dpr data_SENSE3_MB3_noisemap.nii.gz sigma_nmaps.nii.gz N_nmaps.nii.gz mask_nmaps.nii.gz --noise_maps -f --fast_median -m maxlk",
    "dpr dwi_1_8.nii.gz sigma.nii.gz N.nii.gz mask.nii.gz -v",
    "dpr dwi_1_8.nii.gz sigma.nii.gz N.nii.gz mask.nii.gz -m maxlk -f --ncores 4",
    "dpr dwi_1_8.nii.gz sigma_maxlk.nii.gz N_maxlk.nii.gz mask_maxlk.nii.gz -m maxlk --size 3 -f -v --axis 0",
]


@pytest.mark.parametrize("command", commands_dpr)
def test_script_dpr(command):
    subprocess.run([command], shell=True, cwd=cwd, check=True)


@pytest.mark.parametrize("command", commands_dpr_graph)
def test_script_dpr_graph(command):
    subprocess.run([command], shell=True, cwd=cwd, check=True)

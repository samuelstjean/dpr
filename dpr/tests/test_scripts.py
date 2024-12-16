import subprocess
import pytest

from pathlib import Path

cwd = Path(__file__).parents[2] / Path("datasets")
commands_dpr = [
    ('dpr', 'af_left_AFD.txt', 'af_left_AFD_realigned.txt', '--exploredti', '--do_graph', '-f', '-v', '--points', '75'),
    ('dpr', 'af_left_AFD.txt', 'af_left_AFD_realigned.csv', '--exploredti', '-f')
]

commands_dpr_graph = [
    ('dpr_make_fancy_graph', 'af_left_pval_unaligned.txt', 'af_left_coordinates.txt', 'af_left_truncated_coordinates.txt', 'af_left_average_coordinates.txt', '0,2', 'pvals_unaligned.png', '--title', "p-values before realignment", '-f'),
    ('dpr_make_fancy_graph', 'af_left_pval_realigned.txt', 'af_left_coordinates.txt', 'af_left_truncated_coordinates.txt', 'af_left_average_coordinates.txt', '0,2', 'pvals_realigned.png', '-f', '-v', '--labelx', '"label X"', '--labely', 'Label Y', '--dpi', '100')
]

@pytest.mark.parametrize("command", commands_dpr)
def test_script_dpr(command):
    subprocess.run(command, cwd=cwd, check=True)

@pytest.mark.parametrize("command", commands_dpr_graph)
def test_script_dpr_graph(command):
    subprocess.run(command, cwd=cwd, check=True)

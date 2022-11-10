# -*- coding: utf-8 -*-

# This file was automatically generated by parselglossy on 2022-12-03
# Editing is *STRONGLY DISCOURAGED*

from copy import deepcopy
import json
from pathlib import Path
from typing import Optional, Union

from .plumbing import lexer
from .plumbing.utils import JSONDict, as_complex, path_resolver, copier
from .plumbing.validation import validate_from_dicts


def lex(
    infile: Union[str, Path],
    ir_file: Optional[Union[str, Path]] = None,
) -> JSONDict:
    """Run grammar of choice on input string.

    Parameters
    ----------
    infile : Union[str, Path]
        The string to be parsed.
    ir_file : Optional[Union[str, Path]]
        File to write intermediate representation to (JSON format).
        None by default, which means file is not written out.

    Returns
    -------
    The contents of the input string as a dictionary.
    """

    infile = path_resolver(infile)
    if ir_file is not None:
        ir_file = path_resolver(ir_file)

    with infile.open("r") as f:
        ir = lexer.lex_from_str(in_str=f.read(), ir_file=ir_file)

    return ir


def validate(
    ir_in: Union[str, Path, JSONDict],
    fr_file: Optional[Union[str, Path]] = None,
) -> JSONDict:
    """Validate intermediate representation into final representation.

    Parameters
    ----------
    ir_in : Union[str, Path, JSONDict]
        The file (JSON format) or ``JSONDict`` with the intermediate representation.
    fr_file : Optional[Union[str, Path]]
        File to write final representation to (JSON format).
        None by default, which means file is not written out.

    Returns
    -------
    The validated input as a dictionary.
    """

    if isinstance(ir_in, (str, Path)):
        infile = path_resolver(ir_in)
        with infile.open("r") as f:
            ir = json.load(f, object_hook=as_complex)
    elif isinstance(ir_in, dict):
        ir = deepcopy(ir_in)
    else:
        raise RuntimeError("Unrecognized type for argument ``ir_in``.")

    if fr_file is not None:
        fr_file = path_resolver(fr_file)

    return validate_from_dicts(ir=ir, template=stencil(), fr_file=fr_file)


def parse(
    infile: Union[str, Path],
    outfile: Optional[Union[str, Path]] = None,
    dump_ir: bool = False,
) -> JSONDict:
    """Parse input file.

    Parameters
    ----------
    infile : Union[str, Path]
        The input file to be parsed.
    outfile : Optional[Union[str, Path]]
        The output file.
        Defaults to ``None``, which means writing to ``<infile>_fr.json``.
    dump_ir : bool
        Whether to write out the intermediate representation to file (JSON format).
        False by default. If true the filename if ``<infile>_ir.json``

    Returns
    -------
    The validated input as a dictionary.
    """

    infile = path_resolver(infile)

    ir_file = None  # type: Optional[Path]
    if dump_ir:
        ir_file = path_resolver(Path(infile).stem + "_ir.json")

    ir = lex(infile=infile, ir_file=ir_file)

    if outfile is not None:
        outfile = path_resolver(outfile)

    return validate_from_dicts(ir=ir, template=stencil(), fr_file=outfile)


def stencil() -> JSONDict:
    return {   'keywords': [   {   'name': 'world_prec',
                        'predicates': ['1.0e-10 < value < 1.0'],
                        'type': 'float'},
                    {   'default': -1,
                        'name': 'world_size',
                        'predicates': ['value <= 10'],
                        'type': 'int'},
                    {   'default': 'bohr',
                        'name': 'world_unit',
                        'predicates': ['value.lower() in ["bohr", "angstrom"]'],
                        'type': 'str'},
                    {   'default': [0.0, 0.0, 0.0],
                        'name': 'world_origin',
                        'predicates': ['len(value) == 3'],
                        'type': 'List[float]'}],
    'sections': [   {   'keywords': [   {   'default': -1.0,
                                            'name': 'exchange_prec',
                                            'type': 'float'},
                                        {   'default': -1.0,
                                            'name': 'helmholtz_prec',
                                            'type': 'float'},
                                        {   'default': "user['world_prec']",
                                            'name': 'poisson_prec',
                                            'predicates': [   '1.0e-10 < value '
                                                              '< 1.0'],
                                            'type': 'float'},
                                        {   'default': "user['world_prec']",
                                            'name': 'nuclear_prec',
                                            'predicates': [   '1.0e-10 < value '
                                                              '< 1.0'],
                                            'type': 'float'}],
                        'name': 'Precisions'},
                    {   'keywords': [   {   'default': 0,
                                            'name': 'print_level',
                                            'type': 'int'},
                                        {   'default': False,
                                            'name': 'print_mpi',
                                            'type': 'bool'},
                                        {   'default': 6,
                                            'name': 'print_prec',
                                            'predicates': ['0 < value < 10'],
                                            'type': 'int'},
                                        {   'default': 75,
                                            'name': 'print_width',
                                            'predicates': ['50 < value < 100'],
                                            'type': 'int'},
                                        {   'default': False,
                                            'name': 'print_constants',
                                            'type': 'bool'}],
                        'name': 'Printer'},
                    {   'keywords': [   {   'default': 'plots',
                                            'name': 'path',
                                            'predicates': ["value[-1] != '/'"],
                                            'type': 'str'},
                                        {   'default': 'cube',
                                            'name': 'type',
                                            'predicates': [   'value.lower() '
                                                              "in ['line', "
                                                              "'surf', "
                                                              "'cube']"],
                                            'type': 'str'},
                                        {   'default': [20, 20, 20],
                                            'name': 'points',
                                            'predicates': [   'all(p > 0 for p '
                                                              'in value)',
                                                              'not '
                                                              "(user['Plotter']['type'] "
                                                              "== 'line' and "
                                                              'len(value) < 1)',
                                                              'not '
                                                              "(user['Plotter']['type'] "
                                                              "== 'surf' and "
                                                              'len(value) < 2)',
                                                              'not '
                                                              "(user['Plotter']['type'] "
                                                              "== 'cube' and "
                                                              'len(value) < '
                                                              '3)'],
                                            'type': 'List[int]'},
                                        {   'default': [0.0, 0.0, 0.0],
                                            'name': 'O',
                                            'predicates': ['len(value) == 3'],
                                            'type': 'List[float]'},
                                        {   'default': [1.0, 0.0, 0.0],
                                            'name': 'A',
                                            'predicates': ['len(value) == 3'],
                                            'type': 'List[float]'},
                                        {   'default': [0.0, 1.0, 0.0],
                                            'name': 'B',
                                            'predicates': ['len(value) == 3'],
                                            'type': 'List[float]'},
                                        {   'default': [0.0, 0.0, 1.0],
                                            'name': 'C',
                                            'predicates': ['len(value) == 3'],
                                            'type': 'List[float]'}],
                        'name': 'Plotter'},
                    {   'keywords': [   {   'default': False,
                                            'name': 'numerically_exact',
                                            'type': 'bool'},
                                        {   'default': 10000,
                                            'name': 'shared_memory_size',
                                            'type': 'int'},
                                        {   'default': False,
                                            'name': 'share_nuclear_potential',
                                            'type': 'bool'},
                                        {   'default': False,
                                            'name': 'share_coulomb_potential',
                                            'type': 'bool'},
                                        {   'default': False,
                                            'name': 'share_xc_potential',
                                            'type': 'bool'},
                                        {   'default': -1,
                                            'name': 'bank_size',
                                            'type': 'int'}],
                        'name': 'MPI'},
                    {   'keywords': [   {   'default': -1,
                                            'name': 'order',
                                            'type': 'int'},
                                        {   'default': 'interpolating',
                                            'name': 'type',
                                            'predicates': [   'value.lower() '
                                                              'in '
                                                              "['interpolating', "
                                                              "'legendre']"],
                                            'type': 'str'}],
                        'name': 'Basis'},
                    {   'keywords': [   {   'default': 'abgv_55',
                                            'name': 'kinetic',
                                            'type': 'str'},
                                        {   'default': 'abgv_00',
                                            'name': 'h_b_dip',
                                            'type': 'str'},
                                        {   'default': 'abgv_00',
                                            'name': 'h_m_pso',
                                            'type': 'str'},
                                        {   'default': 'abgv_00',
                                            'name': 'zora',
                                            'type': 'str'}],
                        'name': 'Derivatives'},
                    {   'keywords': [   {   'default': 0,
                                            'name': 'charge',
                                            'type': 'int'},
                                        {   'default': 1,
                                            'name': 'multiplicity',
                                            'predicates': ['value > 0'],
                                            'type': 'int'},
                                        {   'default': False,
                                            'name': 'translate',
                                            'type': 'bool'},
                                        {'name': 'coords', 'type': 'str'}],
                        'name': 'Molecule'},
                    {   'keywords': [   {   'name': 'method',
                                            'predicates': [   'value.lower() '
                                                              "in ['core', "
                                                              "'hartree', "
                                                              "'hf', "
                                                              "'hartreefock', "
                                                              "'hartree-fock', "
                                                              "'dft', 'lda', "
                                                              "'svwn3', "
                                                              "'svwn5', 'pbe', "
                                                              "'pbe0', "
                                                              "'bpw91', "
                                                              "'bp86', "
                                                              "'b3p86', "
                                                              "'b3p86-g', "
                                                              "'blyp', "
                                                              "'b3lyp', "
                                                              "'b3lyp-g', "
                                                              "'olyp', 'kt1', "
                                                              "'kt2', 'kt3']"],
                                            'type': 'str'},
                                        {   'default': True,
                                            'name': 'restricted',
                                            'type': 'bool'},
                                        {   'default': 'none',
                                            'name': 'relativity',
                                            'predicates': [   'value.lower() '
                                                              "in ['none', "
                                                              "'zora', "
                                                              "'nzora']"],
                                            'type': 'str'},
                                        {   'default': 'none',
                                            'name': 'environment',
                                            'predicates': [   'value.lower() '
                                                              "in ['none', "
                                                              "'pcm']"],
                                            'type': 'str'}],
                        'name': 'WaveFunction'},
                    {   'keywords': [   {   'default': True,
                                            'name': 'include_nuclear',
                                            'type': 'bool'},
                                        {   'default': True,
                                            'name': 'include_coulomb',
                                            'type': 'bool'},
                                        {   'default': True,
                                            'name': 'include_xc',
                                            'type': 'bool'}],
                        'name': 'ZORA'},
                    {   'keywords': [   {   'default': 0.0,
                                            'name': 'density_cutoff',
                                            'type': 'float'},
                                        {   'default': ' ',
                                            'name': 'functionals',
                                            'type': 'str'},
                                        {   'default': "not(user['WaveFunction']['restricted'])",
                                            'name': 'spin',
                                            'type': 'bool'}],
                        'name': 'DFT'},
                    {   'keywords': [   {   'default': True,
                                            'name': 'dipole_moment',
                                            'type': 'bool'},
                                        {   'default': False,
                                            'name': 'quadrupole_moment',
                                            'type': 'bool'},
                                        {   'default': False,
                                            'name': 'polarizability',
                                            'type': 'bool'},
                                        {   'default': False,
                                            'name': 'magnetizability',
                                            'type': 'bool'},
                                        {   'default': False,
                                            'name': 'nmr_shielding',
                                            'type': 'bool'},
                                        {   'default': False,
                                            'name': 'geometric_derivative',
                                            'type': 'bool'},
                                        {   'default': False,
                                            'name': 'plot_density',
                                            'type': 'bool'},
                                        {   'default': [],
                                            'name': 'plot_orbitals',
                                            'type': 'List[int]'}],
                        'name': 'Properties'},
                    {   'keywords': [   {   'default': [],
                                            'name': 'electric_field',
                                            'predicates': [   'len(value) == 0 '
                                                              'or len(value) '
                                                              '== 3'],
                                            'type': 'List[float]'}],
                        'name': 'ExternalFields'},
                    {   'keywords': [   {   'default': [0.0],
                                            'name': 'frequency',
                                            'type': 'List[float]'}],
                        'name': 'Polarizability'},
                    {   'keywords': [   {   'default': False,
                                            'name': 'nuclear_specific',
                                            'type': 'bool'},
                                        {   'default': [-1],
                                            'name': 'nucleus_k',
                                            'type': 'List[int]'}],
                        'name': 'NMRShielding'},
                    {   'keywords': [   {   'default': 'initial_guess/mrchem.bas',
                                            'name': 'guess_basis',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/mrchem.mop',
                                            'name': 'guess_gto_p',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/mrchem.moa',
                                            'name': 'guess_gto_a',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/mrchem.mob',
                                            'name': 'guess_gto_b',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/phi_p',
                                            'name': 'guess_phi_p',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/phi_a',
                                            'name': 'guess_phi_a',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/phi_b',
                                            'name': 'guess_phi_b',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/X_p',
                                            'name': 'guess_x_p',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/X_a',
                                            'name': 'guess_x_a',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/X_b',
                                            'name': 'guess_x_b',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/Y_p',
                                            'name': 'guess_y_p',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/Y_a',
                                            'name': 'guess_y_a',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/Y_b',
                                            'name': 'guess_y_b',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/phi_p',
                                            'name': 'guess_cube_p',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/phi_a',
                                            'name': 'guess_cube_a',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/phi_b',
                                            'name': 'guess_cube_b',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/x_p',
                                            'name': 'guess_cube_x_p',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/x_a',
                                            'name': 'guess_cube_x_a',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/x_b',
                                            'name': 'guess_cube_x_b',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/y_p',
                                            'name': 'guess_cube_y_p',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/y_a',
                                            'name': 'guess_cube_y_a',
                                            'type': 'str'},
                                        {   'default': 'initial_guess/y_b',
                                            'name': 'guess_cube_y_b',
                                            'type': 'str'},
                                        {   'default': 'cube_vectors/',
                                            'name': 'cube_vectors',
                                            'type': 'str'}],
                        'name': 'Files'},
                    {   'keywords': [   {   'default': True,
                                            'name': 'run',
                                            'type': 'bool'},
                                        {   'default': 100,
                                            'name': 'max_iter',
                                            'type': 'int'},
                                        {   'default': 5,
                                            'name': 'kain',
                                            'type': 'int'},
                                        {   'default': 0,
                                            'name': 'rotation',
                                            'type': 'int'},
                                        {   'default': False,
                                            'name': 'localize',
                                            'type': 'bool'},
                                        {   'default': -1.0,
                                            'name': 'energy_thrs',
                                            'type': 'float'},
                                        {   'default': 0.001,
                                            'name': 'guess_prec',
                                            'predicates': [   '1.0e-10 < value '
                                                              '< 1.0'],
                                            'type': 'float'},
                                        {   'default': 12.0,
                                            'name': 'guess_screen',
                                            'type': 'float'},
                                        {   'default': -1.0,
                                            'name': 'start_prec',
                                            'type': 'float'},
                                        {   'default': -1.0,
                                            'name': 'final_prec',
                                            'type': 'float'},
                                        {   'default': 'sad_gto',
                                            'name': 'guess_type',
                                            'predicates': [   'value.lower() '
                                                              "in ['mw', "
                                                              "'chk', 'gto', "
                                                              "'core_sz', "
                                                              "'core_dz', "
                                                              "'core_tz', "
                                                              "'core_qz', "
                                                              "'sad_sz', "
                                                              "'sad_dz', "
                                                              "'sad_tz', "
                                                              "'sad_qz', "
                                                              "'sad_gto', "
                                                              "'cube']"],
                                            'type': 'str'},
                                        {   'default': False,
                                            'name': 'write_checkpoint',
                                            'type': 'bool'},
                                        {   'default': 'checkpoint',
                                            'name': 'path_checkpoint',
                                            'predicates': ["value[-1] != '/'"],
                                            'type': 'str'},
                                        {   'default': False,
                                            'name': 'write_orbitals',
                                            'type': 'bool'},
                                        {   'default': 'orbitals',
                                            'name': 'path_orbitals',
                                            'predicates': ["value[-1] != '/'"],
                                            'type': 'str'},
                                        {   'default': '10 * '
                                                       "user['world_prec']",
                                            'name': 'orbital_thrs',
                                            'type': 'float'}],
                        'name': 'SCF'},
                    {   'keywords': [   {   'default': [True, True, True],
                                            'name': 'run',
                                            'type': 'List[bool]'},
                                        {   'default': 100,
                                            'name': 'max_iter',
                                            'type': 'int'},
                                        {   'default': 5,
                                            'name': 'kain',
                                            'type': 'int'},
                                        {   'default': -1.0,
                                            'name': 'property_thrs',
                                            'type': 'float'},
                                        {   'default': -1.0,
                                            'name': 'start_prec',
                                            'type': 'float'},
                                        {   'default': -1.0,
                                            'name': 'final_prec',
                                            'type': 'float'},
                                        {   'default': 0.001,
                                            'name': 'guess_prec',
                                            'predicates': [   '1.0e-10 < value '
                                                              '< 1.0'],
                                            'type': 'float'},
                                        {   'default': 'none',
                                            'name': 'guess_type',
                                            'predicates': [   'value.lower() '
                                                              "in ['none', "
                                                              "'chk', 'mw', "
                                                              "'cube', "
                                                              "'random']"],
                                            'type': 'str'},
                                        {   'default': 0,
                                            'name': 'seed',
                                            'type': 'int'},
                                        {   'default': False,
                                            'name': 'write_checkpoint',
                                            'type': 'bool'},
                                        {   'default': 'checkpoint',
                                            'name': 'path_checkpoint',
                                            'predicates': ["value[-1] != '/'"],
                                            'type': 'str'},
                                        {   'default': False,
                                            'name': 'write_orbitals',
                                            'type': 'bool'},
                                        {   'default': 'orbitals',
                                            'name': 'path_orbitals',
                                            'predicates': ["value[-1] != '/'"],
                                            'type': 'str'},
                                        {   'default': '10 * '
                                                       "user['world_prec']",
                                            'name': 'orbital_thrs',
                                            'type': 'float'},
                                        {   'default': "user['SCF']['localize']",
                                            'name': 'localize',
                                            'type': 'bool'}],
                        'name': 'Response'},
                    {   'name': 'PCM',
                        'sections': [   {   'keywords': [   {   'default': 100,
                                                                'name': 'max_iter',
                                                                'type': 'int'},
                                                            {   'default': True,
                                                                'name': 'dynamic_thrs',
                                                                'type': 'bool'},
                                                            {   'default': 'potential',
                                                                'name': 'optimizer',
                                                                'predicates': [   'value.lower() '
                                                                                  'in '
                                                                                  "['density', "
                                                                                  "'potential']"],
                                                                'type': 'str'},
                                                            {   'default': 'total',
                                                                'name': 'density_type',
                                                                'predicates': [   'value.lower() '
                                                                                  'in '
                                                                                  "['total', "
                                                                                  "'nuclear', "
                                                                                  "'electronic']"],
                                                                'type': 'str'},
                                                            {   'default': "user['SCF']['kain']",
                                                                'name': 'kain',
                                                                'type': 'int'}],
                                            'name': 'SCRF'},
                                        {   'keywords': [   {   'default': 'atoms',
                                                                'name': 'mode',
                                                                'predicates': [   'value.lower() '
                                                                                  'in '
                                                                                  "['atoms', "
                                                                                  "'explicit']"],
                                                                'type': 'str'},
                                                            {   'default': '',
                                                                'name': 'spheres',
                                                                'type': 'str'},
                                                            {   'default': 1.1,
                                                                'name': 'alpha',
                                                                'type': 'float'},
                                                            {   'default': 0.5,
                                                                'name': 'beta',
                                                                'type': 'float'},
                                                            {   'default': 0.2,
                                                                'name': 'sigma',
                                                                'type': 'float'}],
                                            'name': 'Cavity'},
                                        {   'keywords': [   {   'default': 1.0,
                                                                'name': 'epsilon_in',
                                                                'type': 'float'},
                                                            {   'default': 1.0,
                                                                'name': 'epsilon_out',
                                                                'type': 'float'},
                                                            {   'default': 'exponential',
                                                                'name': 'formulation',
                                                                'predicates': [   'value.lower() '
                                                                                  'in '
                                                                                  "['exponential']"],
                                                                'type': 'str'}],
                                            'name': 'Permittivity'}]},
                    {   'keywords': [   {   'default': 78.9451185,
                                            'name': 'hartree2simagnetizability',
                                            'type': 'float'},
                                        {   'default': 137.035999084,
                                            'name': 'light_speed',
                                            'type': 'float'},
                                        {   'default': 1.8897261246257702,
                                            'name': 'angstrom2bohrs',
                                            'type': 'float'},
                                        {   'default': 2625.4996394798254,
                                            'name': 'hartree2kjmol',
                                            'type': 'float'},
                                        {   'default': 627.5094740630558,
                                            'name': 'hartree2kcalmol',
                                            'type': 'float'},
                                        {   'default': 27.211386245988,
                                            'name': 'hartree2ev',
                                            'type': 'float'},
                                        {   'default': 219474.6313632,
                                            'name': 'hartree2wavenumbers',
                                            'type': 'float'},
                                        {   'default': 0.0072973525693,
                                            'name': 'fine_structure_constant',
                                            'type': 'float'},
                                        {   'default': -2.00231930436256,
                                            'name': 'electron_g_factor',
                                            'type': 'float'},
                                        {   'default': 2.5417464739297717,
                                            'name': 'dipmom_au2debye',
                                            'type': 'float'}],
                        'name': 'Constants'}]}

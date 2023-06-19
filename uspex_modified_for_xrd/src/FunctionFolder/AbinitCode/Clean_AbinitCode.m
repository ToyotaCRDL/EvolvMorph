function Clean_AbinitCode(code)
if code == 1
[a,b]= unix('mv CONTCAR CONTCAR_old');
[a,b]= unix('mv OUTCAR OUTCAR_old');
[a,b]= unix('mv OSZICAR OSZICAR_old');
[a,b]= unix('mv POSCAR POSCAR_old');
[a,b]= unix('mv vasprun.xml vasprun.xml_old');   # added on 210414
elseif code ==2
[a,b]=unix('rm *-pbe *.ion* siesta.* INPUT_TMP.* *.psdump *.confpot');
elseif code ==3
[a,b]=unix('rm input optimized.structure');
[a,b]=unix('mv output gulp-old.output');
elseif code ==6
[a,b]=unix('rm fort.12 fort.13 dma_output mol.dmain');
elseif code ==7
[a,b]=unix('mv USPEX-1.cell    old_USPEX-1.cell');
[a,b]=unix('mv USPEX-pos-1.xyz old_USPEX-pos-1.xyz');
[a,b]=unix('rm -f USPEX*  *.uspex  cp2k.inp  lattice.stress_tensor');
elseif code ==8
[a,b]=unix('mv output qE_old.output');
elseif code ==9
[a,b]=unix('mv FHI_output FHI_output_old');
[a,b]=unix('rm relaxation_restart_file.FHIaims');
elseif code == 12
[a,b]=unix('rm -f output *.current_stage *.int* *.key *.make0 *.make *.seq *.tmp *.xyz* *.angles* *.pdb* POSCAR*');
end

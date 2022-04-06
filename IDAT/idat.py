with open(f'{swarm_scripts_dir}/idat_to_ped.swarm', 'w') as f:
    for idat_dir in [p1_idat_dir, p23_idat_dir, p1216_idat_dir, mark_caroline, p411re_idat_dir, p17_idat_dir]:
        idat_to_ped_cmd = f'\
{iaap} gencall \
{bpm} \
{egt} \
{ped_dir}/ \
-f {idat_dir}/ \
-p \
-t 16'
        
        f.write(f'{idat_to_ped_cmd}\n')
f.close()


ped_filenames = [p.split('.')[0] for p in glob.glob(f'{basedir}/ped/*.ped')]
map_file = f'{ped_dir}/NeuroBooster_20042459_A2.map'
for filename in ped_filenames:
    ped = f'{filename}.ped'
    out_map = f'{filename}.map'
    if os.path.isfile(ped):
        shutil.copyfile(src=map_file, dst=out_map)
    else:
        print(f'{ped} does not exist!')
        print(f'{out_map} creation cancelled')
        
with open(f'{swarm_scripts_dir}/make_bed.swarm', 'w') as f:
    for filename in ped_filenames:
        ped = f'{filename}'
        make_bed_cmd = f'\
plink \
--file {ped} \
--make-bed \
--out {filename}'

        f.write(f'{make_bed_cmd}\n')
f.close()


with open(f"{bed_dir}/merge_bed.list", 'w') as f:
    for filename in ped_filenames:
        bed = f'{filename}'
        f.write(f'{bed}\n')
f.close()

with open(f"{swarm_scripts_dir}/merge.swarm", 'w') as f:

    plink_merge_cmd = f'\
plink \
--merge-list {bed_dir}/merge_bed.list \
--make-bed \
--out {bed_dir}/ibx'
    f.write(f"{plink_merge_cmd}")
f.close()

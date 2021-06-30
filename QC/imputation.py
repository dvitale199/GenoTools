import subprocess
import argparse
import os
import requests
import json
import time
import glob
import shutil

# local imports
from QC.utils import shell_do

def impute_data_prep(geno_path, out_path, ref_panel, check_bim_pl):

    '''
    info here:
    https://www.well.ox.ac.uk/~wrayner/tools/
    '''
    
    out_dir = out_path.rsplit('/', 1)[0]
    workdir = os.getcwd()
    os.chdir(out_dir)
    
    plink_files = [f'{geno_path}.{suffix}' for suffix in ['bed','bim','fam']]
    ref_files = [ref_panel, check_bim_pl]
    
    cp_files = plink_files + ref_files
    
    for file in cp_files:
        if os.path.isfile(file):
            shutil.copy(file, out_dir)
        else:
            print(f'ERROR: {file} does not exist!')
            
    
    geno_path2 = geno_path.split('/')[-1]
    out_path2 = out_path.split('/')[-1]
    ref_panel2 = ref_panel.split('/')[-1]
    check_bim_pl2 = check_bim_pl.split('/')[-1]
    
    
    plink1 = f'plink --bfile {geno_path2} --freq --out {out_path2}'
    check_bim_cmd = f'perl {check_bim_pl2} -b {geno_path2}.bim -f {out_path2}.frq -r {ref_panel2} -h'
    bash1 = 'sh Run-plink.sh'

    cmds = [plink1, check_bim_cmd, bash1]

    for cmd in cmds:

        shell_do(cmd, log=True, return_log=True)


    
    mk_vcf_cmds = [f'plink --bfile {geno_path2}-updated-chr{str(i)} --recode vcf --chr {str(i)} --out {out_path2}_chr{str(i)}' for i in range(1,24)]   

    for cmd in mk_vcf_cmds:

        shell_do(cmd)

    # then sort and zip
    sort_zip_cmds = [f'vcf-sort {out_path2}_chr{str(i)}.vcf | bgzip -c > {out_path2}_pre_impute_chr{str(i)}.vcf.gz' for i in range(1,24)]

    for cmd in sort_zip_cmds:

        subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
    
    os.chdir(workdir)
    
    vcf_outpaths = [f'{out_path}_pre_impute_chr{str(i)}.vcf.gz' for i in range(1,24)]
    return {'vcfs': vcf_outpaths}



def check_impute_status(token, job_id):
        
    # imputation server url
    url = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'

    # add token to header (see authentication)
    headers = {'X-Auth-Token' : token }

    # get all jobs
    r = requests.get(url + "/jobs", headers=headers)
    if r.status_code != 200:
        raise Exception('GET /jobs/ {}'.format(r.status_code))

    status = r.json()
    for stat in status['data']:
        if stat['id'] == job_id:
            if stat['state'] == 1:
                print("Launching Job:", stat['id'])
            elif stat['state'] == 2:
                print("Running Job:", stat['id'])
            elif stat['state'] == 3:
                print(stat['id'], "returned state '3', have a look at jobs on the web front for more information")
            elif stat['state'] == 5:
                print(stat['id'], "has failed. consult docs on data input to ensure your vcfs are correct")
            elif stat['state'] == 4:
                print(stat['id'], "COMPLETED!")

            return stat['state']

        else:
            pass

        
def pull_imputed_data(out_path, token, job_id, password):
    
    workdir = os.getcwd()
    os.chdir(out_path)
    
    # imputation server url
    url = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'

    # add token to header (see authentication)
    headers = {'X-Auth-Token' : token }

    r = requests.get(f'{url}/jobs/{job_id}', headers=headers)
    if r.status_code != 200:
        raise Exception(f'GET /jobs/ {r.status_code}')

    output_json = r.json()

    hashes_dict = {output_json['outputParams'][i]['id'] : output_json['outputParams'][i]['hash'] for i in range(len(output_json['outputParams']))}

    # run a curl for each
    curls = [f'curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/{str(key)}/{str(hashes_dict[key])} | bash' for key in hashes_dict]
    
    for curl in curls:
        print(f"Curling output data with the following command: {curl}")
        subprocess.run(curl, shell=True)
    print() 
    print("Finished Pulling Imputed Data!")
    print()
    
    os.chdir(workdir)
    

def submit_job(vcf_list,
               password, token=None):
    # test topmed server
    url = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'
    
    # add token to header (see Authentication)
    headers = {'X-Auth-Token' : token}

    open_vcfs = [open(vcf, 'rb') for vcf in vcf_list]

    files = set([('input-files-upload', vcf) for vcf in open_vcfs])
    # files = {'input-files': open(vcf_list[0], 'rb')}

    data = {'input-mode' : 'imputation',
            'input-files-source': 'file-upload',
            'input-password': password,
            'input-refpanel': 'apps@topmed-r2@1.0.0',
            'input-phasing': 'eagle',
            'input-population': 'all'}

    r = requests.post(url + "/jobs/submit/imputationserver", files=files, headers=headers, data=data)
    if r.status_code != 200:
        raise Exception('POST /jobs/submit/imputationserver {}'.format(r.status_code))

    job_id = r.json()['id']
    message = r.json()['message']
    print(job_id, message)
    
    print('***************************')
    print('* * * * * * * * * * * * * *')
    
    return r.json()




def run_auto_imputation(vcf_list, out_path, token, password='imputer'):
    

#     impute_data = impute_data_prep(geno_path, temp_path, ref_panel, check_bim_pl)    
#     job_json = submit_job(impute_data['vcfs'], password=password, token=token)

    job_json = submit_job(vcf_list, password=password, token=token)
    
    job_id = job_json['id']
    
    imp_state = 0
    while imp_state < 3:
        time.sleep(600)
        os.system('clear')
        imp_state = check_impute_status(token, job_id)

        if imp_state == 4:
            print("Pulling Completed Data from Imputation Server!")
            pull_imputed_data(out_path=out_path, token=token, job_id=job_id, password=password)
    
    return job_json

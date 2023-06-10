import os

# log_files = ['/Users/kuznetsovn2/Desktop/GitHub/GenoTools/data/GP2_QC_round3_S4_umap_linearsvc_adjusted_labels.txt', 
#             '/Users/kuznetsovn2/Desktop/GitHub/GenoTools/data/GP2_QC_round3_S4_umap_linearsvc_predicted_labels.txt']

# with open(f'tester.log', "a+") as new_created_file:
#     for name in log_files:
#         steps = f'hello {name}'
#         with open(name) as file:
#             copy_over = file.readlines()
#             # for line in file:
#             new_created_file.seek(0)
#             new_created_file.write(steps)
#             new_created_file.writelines(copy_over)
#             new_created_file.write("\n")

#     tests = new_created_file.readlines()
#     new_created_file.truncate()
#     print(tests)


# def clean_log(concat_log):
#     # to get step name: can subtract bfile from out file
#     # for line in concat_log: # it works! get line: "Step:", keep track of each log with "step"
#     #     print(line)
#     #     break

#     error_index = []
#     warning_index = []
#     for line_num in range(len(concat_log)):
#         if "error" in concat_log[line_num].lower():
#             error_index.append(line_num)
#         if "warning" in concat_log[line_num].lower():
#             warning_index.append(line_num)
#     print(error_index)
#         # print(concat_log[line_num])

def replace_all(text, dict):
    for i, j in dict.items():
        text = text.replace(i, j)
    return text

def process_log(out_dir, mega_log):
    exclude = ['Hostname', 'Working directory', 'Intel', 'Start time', 'Random number seed', 'RAM detected', 'threads', 'thread', 
    'written to', 'done.', 'End time:', 'Writing', '.bed', '.bim', '.fam', '.id', '.hh', '.sexcheck', '.psam', '-bit',
    '.pvar', '.pgen', '.in', '.out', '.het', '.missing', '.snplist', '.kin0', '.eigenvec', '.eigenval', '(--maf/', 'Step:']

    step_indices = [i for i, s in enumerate(mega_log) if 'Step:' in s]
    process_indices = [i for i, s in enumerate(mega_log) if 'Process:' in s]
    bfile_indices = [i for i, s in enumerate(mega_log) if '--bfile' in s]
    pfile_indices = [i for i, s in enumerate(mega_log) if '--pfile' in s]

    bfile_indices.extend(pfile_indices)
    bfile_indices.sort()

    step_indices.append(len(mega_log))
    start = 0
    stop = 1

    with open("data/cleaned_log.gt", "w") as f:
        while start < len(step_indices)-1:
            # step names
            step_line = mega_log[step_indices[start]]
            process_line = mega_log[process_indices[start]]
            process_name = process_line.split('_')[0].replace('Process: ', '')
            bfile_line = mega_log[bfile_indices[start]]
            bfile_name = bfile_line.split()[1].replace(out_dir, "")

            if process_name == 'plink':
                process_name = 'pca'

            step_name = step_line.split(f'{process_name}')[-1]
            step_name = (process_name + step_name).replace('.log', '')

            if bfile_name.split("_")[-1].isupper():
                ancestry_check = bfile_name.split("_")[-1]

            f.write(f'Step: {step_name}')
            f.write(f'Ancestry: {ancestry_check}\n')
            
            for i in range(step_indices[start], step_indices[stop]):
                if not any([x in mega_log[i].strip('\n') for x in exclude]):
                    f.write(mega_log[i])
            start += 1
            stop += 1

# def process_log(out_dir, mega_log):
#     exclude = ['Hostname', 'Working directory', 'Start time', 'Random number seed', 'RAM detected', 'threads', 'thread', 
#     'written to', 'done.', 'End time:', 'Writing', '.bed', '.bim', '.fam', '.id', '.hh', '.sexcheck', '.psam',
#     '.pvar', '.pgen', '.in', '.out', '.het', '.missing', '.snplist', '.kin0', '.eigenvec', '.eigenval', '(--maf/', 'Step:']

#     step_indices = [i for i, s in enumerate(mega_log) if 'Step:' in s]
#     bfile_indices = [i for i, s in enumerate(mega_log) if '--bfile' in s]
#     pfile_indices = [i for i, s in enumerate(mega_log) if '--pfile' in s]

#     bfile_indices.extend(pfile_indices)
#     bfile_indices.sort()
#     step_indices.append(len(mega_log))
#     start = 0
#     stop = 1

#     with open("data/cleaned_log.gt", "w") as f:
#         while start < len(step_indices)-1:
#             # step names
#             step_line = mega_log[step_indices[start]]
#             bfile_line = mega_log[bfile_indices[start]]
#             bfile_name = bfile_line.split()[1].replace(out_dir, "")
#             step_name = step_line.split()[1]

#             replace_dict = {out_dir: "", bfile_name: "", "_": " ", ".logPLINK": ""}
#             step = replace_all(step_name, replace_dict)

#             if bfile_name.split("_")[-1].isupper():
#                 ancestry_check = bfile_name.split("_")[-1]
            
#             if not step: # when bfile = out file
#                 f.write(f'Step: {bfile_name.replace("_", " ")} \n')
#             else:
#                 f.write(f'Step: {ancestry_check} {step} \n')
            
#             for i in range(step_indices[start], step_indices[stop]):
#                 if not any([x in mega_log[i].strip('\n') for x in exclude]):
#                     f.write(mega_log[i])
#             start += 1
#             stop += 1
        

# def process_log(out_path, mega_log):
#     exclude = ['Hostname', 'Working directory', 'Start time', 'Random number seed', 'RAM detected', 'threads', 'thread', 
#     'written to', 'done.', 'End time:', 'Writing', '.bed', '.bim', '.fam', '.id', '.hh', '.sexcheck', '.psam',
#     '.pvar', '.pgen', '.in', '.out', '.het', '.missing', '.snplist', '.kin0', '.eigenvec', '.eigenval', '(--maf/']

#     with open("data/cleaned_log.gt", "w") as f:
#         for line in mega_log:
#             if not any([x in line.strip('\n') for x in exclude]):
#                 f.write(line)

with open('data/all_plink_logs_new.gtlog', 'r') as file:
    # clean_log(file.readlines())
    # out_dir = os.path.dirname(os.path.abspath('data/all_plink_logs_SHORT.gtlog')) # simulate obtaining processing directory
    # print(out_dir)

    out_dir = '/data/CARD_AA/users/kuznetsovn2/'
    process_log(out_dir, file.readlines()) # separates each line
    # we can truncate here so can rewrite into same name - in real thing: send file name or outpath to next function

# print(log_files[0].split('_')[-1])

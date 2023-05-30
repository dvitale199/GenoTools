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


def clean_log(concat_log):
    # to get step name: can subtract bfile from out file
    # for line in concat_log: # it works! get line: "Step:", keep track of each log with "step"
    #     print(line)
    #     break

    error_index = []
    warning_index = []
    for line_num in range(len(concat_log)):
        if "error" in concat_log[line_num].lower():
            error_index.append(line_num)
        if "warning" in concat_log[line_num].lower():
            warning_index.append(line_num)
    print(error_index)
        # print(concat_log[line_num])

with open('data/test_genotools_all_plink_logs.txt', 'r') as file:
    clean_log(file.readlines())
    # we can truncate here so can rewrite into same name - in real thing: send file name or outpath to next function

# print(log_files[0].split('_')[-1])

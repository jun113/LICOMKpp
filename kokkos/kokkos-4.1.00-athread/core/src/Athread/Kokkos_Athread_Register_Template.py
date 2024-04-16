import os
import sys

def athread_search_user_functor (search_dir):
	list_functor_path = []
	for root, dirs, files in os.walk(search_dir):
		for file in files:
			if(file[-4:] == '''.hpp'''):
				with open(os.path.join(root, file), 'r') as f:
					user_hpp_code = f.read()
					if (user_hpp_code.find("KOKKOS_REGISTER_FOR_") != -1):
						list_functor_path.append(os.path.join(root, file))
					if (user_hpp_code.find("KOKKOS_REGISTER_REDUCE_") != -1):
						list_functor_path.append(os.path.join(root, file))
	return list_functor_path

def athread_read_user_functor (list_functor_path):
	list_register_info = []

	for functor_path in list_functor_path:
		with open(functor_path, 'r') as f:
			for line in f.readlines():
				index_RF = line.find('KOKKOS_REGISTER_FOR_')
				# TODO multiline comment and macros condition
				if (index_RF != -1):
					index_slash = line.find('''//''')
					if ((index_slash == -1) or \
						 ((index_slash != -1) and (index_slash > index_RF))):
 
						index_param_l = line.find('(')
						index_param_r = line.find(')')
						index_param_c = line.find(',', index_param_l)
 
						function = line[index_param_l + 1 : index_param_c].strip()
						functor  = line[index_param_c + 1 : index_param_r].strip()
 
						index_dim = line[0 : index_param_l].rfind('D')
						src_dim = "_" + line[index_dim - 1 : index_dim + 1]
						functor_info = {
							"key": functor + src_dim[1],
							"function": function + src_dim,
							"dim": int(src_dim[1])
						}
						list_register_info.append(functor_info)
				index_RF = line.find('KOKKOS_REGISTER_REDUCE_')
				# TODO multiline comment and macros condition
				if (index_RF != -1):
					index_slash = line.find('''//''')
					if ((index_slash == -1) or \
						 ((index_slash != -1) and (index_slash > index_RF))):
 
						index_param_l = line.find('(')
						index_param_r = line.find(')')
						index_param_c = line.find(',', index_param_l)
 
						function = line[index_param_l + 1 : index_param_c].strip()
						functor  = line[index_param_c + 1 : index_param_r].strip()
 
						index_dim = line[0 : index_param_l].rfind('D')
						src_dim = "_" + line[index_dim - 1 : index_dim + 1]
						functor_info = {
							"key": functor + src_dim[1],
							"function": function + src_dim,
							"dim": int(src_dim[1])
						}
						list_register_info.append(functor_info)

	return list_register_info

def athread_write_kernel_launch_slave_source_code ( \
				target_file, list_reg_info, list_functor_path):

	# include head code
	list_inc_code = []
	for functor_path in list_functor_path:
		list_inc_code.append('''#include "''' + functor_path + '''"\n''')

	# register code
	# list_reg_code = ["  using node = AthradRegisterFunctionListNode;\n", \
	# 											"  reg_func_list_node = new node;\n",
	# 											"  node* curr_node = reg_func_list_node;\n\n"
	# 								]
	# list_reg_code = ["  using node = AthradRegisterFunctionListNode;\n", \
	# 											"  reg_func_list_node = (node*)malloc(sizeof(node));\n",
	# 											"  node* curr_node = reg_func_list_node;\n\n"
	# 								]
	list_reg_code = ["  using node = AthradRegisterFunctionListNode;\n", \
												"  reg_func_list_node = (node*)libc_aligned_malloc(sizeof(node));\n",
												"  node* curr_node __attribute__ ((aligned(64))) = reg_func_list_node;\n\n"
									]

	for i in range(len(list_reg_info)):
		register_info = list_reg_info[i]
		# list_reg_code.append("  curr_node->key  = const_cast<char *>(\"" + register_info["key"] + "\");\n")
		list_reg_code.append("  curr_node->key  = arr_char_to_arr_int (\n      \"" + register_info["key"] + "\", curr_node->num_intv16);\n")
											 
											#  const_cast<char *>(\"" + register_info["key"] + "\");\n")
		# list_reg_code.append("  curr_node->key  = \"" + register_info["key"] + "\";\n")
		list_reg_code.append("  curr_node->fp   = " + register_info["function"] + ";\n")
		# list_reg_code.append("  curr_node->dim  = " + str(register_info["dim"]) + ";\n")
		if ((i+1) == len(list_reg_info)):
			list_reg_code.append("  curr_node->next = nullptr;\n")
		else:
			# list_reg_code.append("  curr_node->next = new node;\n\n")
			# list_reg_code.append("  curr_node->next = (node*)malloc(sizeof(node));\n\n")
			list_reg_code.append("  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));\n\n")
			list_reg_code.append("  curr_node       = curr_node->next;\n")

	with open(target_file, 'r') as f:
		old_kernel_launch_slave_code = f.readlines()
		index_inc_start = 0
		index_inc_end   = 0
		index_reg_start = 0
		index_reg_end   = 0

		for index_line in range(len(old_kernel_launch_slave_code)):

			if (old_kernel_launch_slave_code[index_line].find("__INCLUDE_FUNCTOR_HPP_START__") != -1):
				index_inc_start = index_line
				continue
			if (old_kernel_launch_slave_code[index_line].find("__INCLUDE_FUNCTOR_HPP_END__") != -1):
				index_inc_end = index_line
				continue

			if (old_kernel_launch_slave_code[index_line].find("__REGISTER_START_") != -1):
				index_reg_start = index_line
				continue
			if (old_kernel_launch_slave_code[index_line].find("__REGISTER_END_") != -1):
				index_reg_end = index_line
				break
		
		new_kernel_launch_slave_code = \
				old_kernel_launch_slave_code[ : index_inc_start + 1] + \
				list_inc_code + \
				old_kernel_launch_slave_code[index_inc_end : index_reg_start + 1] + \
				list_reg_code + \
				old_kernel_launch_slave_code[index_reg_end: ]

	# with open("1.txt", 'w') as f:

	with open(target_file, 'w') as f:
		f.writelines(new_kernel_launch_slave_code)
	return

if __name__ == "__main__":
	top_dir = sys.argv[1]
	kernel_launch_slave_file = sys.argv[2]

	# print("top dir", top_dir)
	# print("kernel file", kernel_launch_slave_file)

	list_functor_path = athread_search_user_functor(top_dir)
	list_reg_info = athread_read_user_functor(list_functor_path)
	athread_write_kernel_launch_slave_source_code(kernel_launch_slave_file, \
																			list_reg_info, list_functor_path)
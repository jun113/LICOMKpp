#include <cctype>

int str_trim_cmp(char const *cmp_str, char const *tar_str) {
	// INPUT:
	// cmp_str: the pointer of char that may contain spaces
	// tar_str: the pointer of targeted char 
	// RETURN:
	// 0: equation
	// 1: no equation
	// Notes:
	// Like all other functions from <cctype>, the behavior of std::isspace is undefined if the argument's value is neither representable as unsigned char nor equal to EOF. To use these functions safely with plain chars (or signed chars), the argument should first be converted to unsigned char

	if (!(cmp_str == nullptr && tar_str == nullptr) && 
			(cmp_str == nullptr || tar_str == nullptr)) {
		return 1;
	}
  // start and end index of target string
	int t_i = 0;
	while (std::isspace(static_cast<unsigned char>(tar_str[t_i]))) {
  	t_i += 1;
  }
  // index of compared string
  int c_i = 0;
  while (std::isspace(static_cast<unsigned char>(cmp_str[c_i]))) {
    c_i += 1;
  }
  while (tar_str[t_i] != '\0') {
  	if (tar_str[t_i] == cmp_str[c_i]) {
      t_i += 1;
      c_i += 1;
    } else {
      return 1;
    } 
  }
  if (cmp_str[c_i] != '\0' && !std::isspace(static_cast<unsigned char>(cmp_str[c_i]))) {
    return 1;
  }
  return 0;
}
#include "main.h"

//RNG and distribution
std::uniform_int_distribution<uint64_t> uint_dist;
std::random_device rd;
std::mt19937_64 mt(static_cast<uint64_t>(rd()));

//subset sum parameters
uint64_t* ss_array;
uint64_t target;
uint64_t target_outer;
uint64_t e;

//Misc
uint64_t bitmask_n = 1ULL << n;
uint64_t modul_squared;
uint64_t**binom;
uint64_t calls = 0;

//function parameters
uint64_t flavour31, flavour41, flavour11, flavour21;
uint64_t random_offset = 1;
uint64_t func_switch[3];
uint64_t xor_switch;
uint64_t int_val;
uint64_t int_val_rev;
uint64_t inner_coll[4];
uint64_t outer_coll[2];

uint64_t choose(uint64_t n, uint64_t k)
{
	if (k == 0)
		return 1;
	return (n * choose(n - 1, k - 1)) / k;
}

//projection from Z_choose(n,wt) to bitstrings of length n with wt bits set
uint64_t projection(uint64_t a)
{
	int wn = n;
	int wk = wt;
	uint64_t v = 0;
	int set = 0;
	while (wn != 0)
	{
		if (set == wt)
			break;
		else if (wn + set == wt)
		{
			v += (1ULL << (wn - 1));
			wn -= 1;
			set += 1;
		}
		else if (a < binom[wn - 1][wk])
			wn -= 1;
		else
		{
			a -= binom[wn - 1][wk];
			v += (1ULL << (wn - 1));
			wn -= 1;
			wk -= 1;
			set += 1;
		}
	}

	return v;
}

//subset sum function modulo modul, if subtract =1 output R-resulst
uint64_t ss_function(uint64_t x, uint64_t R, bool subtract)
{
	uint64_t bitmask = 1ULL;
	uint64_t result = 0;

	for (int i = 0; i < n; ++i)
	{
		if ((bitmask & x) > 0)
			result += ss_array[i];
		bitmask <<= 1;
	}
	result = result %modul;
	if (subtract)
	{
		if (R >= result)
			return (R - result);
		else
			return (modul - result + R);
	}
	else
		return result;
}

uint64_t ss_function_non_modul(uint64_t x, uint64_t R, bool subtract)
{
	uint64_t bitmask = 1ULL;
	uint64_t result = 0;

	for (int i = 0; i < n; ++i)
	{
		if ((bitmask & x) > 0)
			result += ss_array[i];
		bitmask <<= 1;
	}
	result = result%modul_squared;
	if (subtract)
	{
		if (R >= result)
			return (R - result);
		else
			return (modul_squared - result + R);
	}
	else
		return result;
}

// computes <x,y> mod 2
bool scalar(uint64_t x, uint64_t y)
{
	uint64_t z = x&y;
	int wt = 0;
	while (z)
	{
		z &= (z - 1);
		wt += 1;
	}
	return (wt & 1);
}
// computes <x,y> mod 4
int scalar4(uint64_t x, uint64_t y)
{
	uint64_t z = x&y;
	int wt = 0;
	while (z)
	{
		z &= (z - 1);
		wt += 1;
	}
	return wt % 4;
}

uint64_t inner_ss_function(uint64_t x, uint64_t flavour1, uint64_t flavour2, int func_num)
{
	calls++;
	uint64_t x_prime = (x*flavour1 + flavour2) % modul_p;
	bool b = scalar(x_prime, func_switch[func_num]);
	return ss_function(projection(x_prime), func_num == 0 ? int_val : int_val_rev, b);
}

void change_func_switchs()
{
	func_switch[0] = uint_dist(mt) % bitmask_n;
	func_switch[1] = uint_dist(mt) % bitmask_n;
	func_switch[2] = uint_dist(mt) % bitmask_n;
	xor_switch = uint_dist(mt) % bitmask_n;
}

bool collision_finding_inner(uint64_t sp, uint64_t flavour1, uint64_t flavour2, int func_num)
{
	uint64_t flavour1sp = (flavour1*sp + flavour2) % modul_p;
	uint64_t flavour2sp = (flavour21*sp + flavour11) % modul_p;
	uint64_t slow = 1ULL << 63;
	uint64_t fast = inner_ss_function(sp, flavour1sp, flavour2sp, func_num);
	uint64_t counter = 0, lam = 1, power = 1;

	//first phase: find colliding output
	while (slow != fast) {
		if (lam == power)
		{
			slow = fast;
			power <<= 1;
			lam = 0;
		}
		fast = inner_ss_function(fast, flavour1sp, flavour2sp, func_num);
		lam += 1;
		counter += 1;
		if (counter > (1ULL << ((int)ceil(n / 4.0) + 3)))
			return collision_finding_inner((sp + random_offset) % modul, flavour1, flavour2, func_num);
	}

	//second phase: find colliding inputs
	slow = fast = sp;
	for (uint64_t i = 0; i<lam; ++i)
		fast = inner_ss_function(fast, flavour1sp, flavour2sp, func_num);
	if (sp == fast)
		return  collision_finding_inner((sp + random_offset) % modul, flavour1, flavour2, func_num);

	uint64_t slow_tmp = 0, fast_tmp;
	counter = 0;
	while (slow != fast) {
		slow_tmp = slow;
		slow = inner_ss_function(slow, flavour1sp, flavour2sp, func_num);
		fast_tmp = fast;
		fast = inner_ss_function(fast, flavour1sp, flavour2sp, func_num);
		counter += 1;
		if (counter > (1ULL << ((int)ceil(n / 4.0) + 3)))
			return collision_finding_inner((sp + random_offset) % modul, flavour1, flavour2, func_num);

	}

	inner_coll[0 + func_num * 2] = (slow_tmp*flavour1sp + flavour2sp) % modul_p;
	inner_coll[1 + func_num * 2] = (fast_tmp*flavour1sp + flavour2sp) % modul_p;
	return true;
}

uint64_t outer_ss_fuction(uint64_t x, uint64_t flavour1, uint64_t flavour2, uint64_t flavour3, uint64_t flavour4)
{
	uint64_t x_prime = (x*flavour3 + flavour4) % modul_p;
	bool b = scalar(x_prime, func_switch[2]);
	collision_finding_inner(x_prime, flavour1, flavour2, b);

	int toggle = scalar4(x_prime, xor_switch);
	uint64_t a = 0;
	uint64_t out = 0;
	//perform one operation of add, xor, or, addition of ss function output, based on input (for consistent vectors all that operations yield the same result)
	if (toggle == 0)
		a = projection(inner_coll[0 + b * 2]) + projection(inner_coll[1 + b * 2]);
	else if (toggle == 1)
		a = projection(inner_coll[0 + b * 2]) ^ projection(inner_coll[1 + b * 2]);
	else if (toggle == 2)
		a = projection(inner_coll[0 + b * 2]) | projection(inner_coll[1 + b * 2]);
	else
		out = (ss_function_non_modul(projection(inner_coll[0 + b * 2]), 0, 0) + ss_function_non_modul(projection(inner_coll[1 + b * 2]), 0, 0));
	if (a != 0)
		out = ss_function_non_modul(a, 0, 0);
	//if inner collision forms an inconsistent vector, remap it randomly
	out = (out + (projection(inner_coll[0 + b * 2]) & projection(inner_coll[1 + b * 2]))
		*(flavour4*(inner_coll[0 + b * 2] * 2 + inner_coll[1 + b * 2]) + flavour3)) % modul_squared;
	uint64_t result = ((out - (out%modul)) / modul) % modul;
	if (b == 1) {
		if (target_outer>result)
			return (target_outer - result);
		return modul + (target_outer - result);
	}
	return result;
}

bool collision_finding_outer(uint64_t sp, uint64_t flavour1, uint64_t flavour2, uint64_t flavour3, uint64_t flavour4)
{
	uint64_t flavour3sp = (sp*flavour3 + flavour4) % modul_p;
	uint64_t flavour4sp = (sp*flavour41 + flavour31) % modul_p;
	uint64_t slow = 1ULL << 63;// outer_ss_fuction(sp, flavour1, flavour2, flavour3sp, flavour4sp);
	uint64_t fast_tmp = sp;
	uint64_t fast = outer_ss_fuction(sp, flavour1, flavour2, flavour3sp, flavour4sp);

	uint64_t counter = 0, lam = 1, power = 1;
	//first phase: find cycle length
	while (slow != fast) {
		if (lam == power)
		{
			slow = fast;
			power <<= 1;
			lam = 0;
		}
		fast = outer_ss_fuction(fast, flavour1, flavour2, flavour3sp, flavour4sp);
		lam += 1;
		counter += 1;
		if (counter > (1ULL << ((int)ceil(n / 4.0) + 3)))
			return  collision_finding_outer((sp + random_offset) % modul, flavour1, flavour2, flavour3, flavour4);

	}
	slow = fast = sp;
	for (uint64_t i = 0; i<lam; ++i)
		fast = outer_ss_fuction(fast, flavour1, flavour2, flavour3sp, flavour4sp);
	if (sp == fast)
		return  collision_finding_outer((sp + random_offset) % modul, flavour1, flavour2, flavour3, flavour4);
	//second phase: find colliding inputs
	uint64_t slow_tmp = 0;
	counter = 0;
	while (slow != fast) {
		slow_tmp = slow;
		slow = outer_ss_fuction(slow, flavour1, flavour2, flavour3sp, flavour4sp);
		fast_tmp = fast;
		fast = outer_ss_fuction(fast, flavour1, flavour2, flavour3sp, flavour4sp);
		counter += 1;
		if (counter >(1ULL << (((int)ceil(n / 4.0)) + 3)))
			return  collision_finding_outer((sp + random_offset) % modul, flavour1, flavour2, flavour3, flavour4);

	}

	outer_coll[0] = (flavour3sp*slow_tmp + flavour4sp) % modul_p;
	outer_coll[1] = (flavour3sp*fast_tmp + flavour4sp) % modul_p;
	return true;
}

bool search_solution()
{
	uint64_t sp, flavour1, flavour2, flavour3, flavour4;
	uint64_t e_prime = 0;
	uint64_t final_value = 0;
	do
	{
		change_func_switchs();
		do
		{
			int_val = uint_dist(mt) % modul;
			int_val_rev = ((target%modul) + modul - int_val) % modul;
		} while (int_val + int_val_rev >= modul);
		sp = uint_dist(mt) % modul;
		flavour1 = uint_dist(mt) % modul_p;
		flavour2 = uint_dist(mt) % modul_p;
		flavour3 = uint_dist(mt) % modul_p;
		flavour4 = uint_dist(mt) % modul_p;
		flavour11 = uint_dist(mt) % modul_p;
		flavour21 = uint_dist(mt) % modul_p;
		flavour31 = uint_dist(mt) % modul_p;
		flavour41 = uint_dist(mt) % modul_p;

		collision_finding_outer(sp, flavour1, flavour2, flavour3, flavour4);
		int scalar0 = scalar(outer_coll[0], func_switch[2]);
		int scalar1 = scalar(outer_coll[1], func_switch[2]);
		if (scalar0 == scalar1)
			continue;

		collision_finding_inner(outer_coll[0], flavour1, flavour2, scalar0);
		collision_finding_inner(outer_coll[1], flavour1, flavour2, scalar1);
		e_prime = projection(inner_coll[0]) | projection(inner_coll[1])
			| projection(inner_coll[2]) | projection(inner_coll[3]);
		final_value = ss_function_non_modul(e_prime, 0, 0);

	} while (final_value != target);
	return true;
}

void init_subset_sum()
{
	for (int i = 0; i < n; ++i)
		ss_array[i] = uint_dist(mt) % modul_squared;

	int w = 0, a;
	e = 0;
	while (w < n / 2) {
		a = uint_dist(mt) % n;
		if ((e&(1ULL << a)) == 0) {
			e ^= (1ULL << a);
			w += 1;
		}
	}
	target = ss_function_non_modul(e, 0, 0) % modul_squared;
	target_outer = ((target - (target%modul)) / modul) % modul;
}
int main()
{
	//precompute n choose k for bijection
	binom = new uint64_t*[n + 1];
	for (int i = 0; i < n + 1; ++i)
	{
		binom[i] = new uint64_t[wt + 1];
		for (int j = 0; j<wt + 1; ++j)
			binom[i][j] = choose(i, j);
	}
	modul = 2017;// binom[n][wt];
	ss_array = new uint64_t[n];
	modul_squared = (uint64_t)pow(modul, 2);

	//initialize subset sum instance
	init_subset_sum();

	//test performance
	int iterations = 30;

	uint64_t total = 0;
	for (int i = 0; i < iterations; ++i)
	{
		init_subset_sum();
		search_solution();

		std::cout << "\n" << i + 1 << ". done with " << log2(calls);
		total += calls;

		calls = 0;
		fflush(stdout);
	}
	std::cout << "\navg calls= " << log2(total*1.0 / iterations)<<"\n";
	
}
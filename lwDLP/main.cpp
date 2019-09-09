#include "main.h"
//RNG and distribution
std::uniform_int_distribution<uint64_t> uint_dist;
std::random_device rd;
std::mt19937_64 mt(static_cast<uint64_t>(rd()));
std::mt19937_64 prg(static_cast<uint64_t>(rd()));


//Misc
int n =1;
uint64_t bitmask_n = 1ULL << n;
uint64_t**binom;
uint64_t calls = 0;
uint64_t prob;
int max_digits = 10;

//function parameters
uint64_t flavour1, flavour2;
uint64_t func_switch;

//array for collision
uint64_t coll[2];

//dlog instance parameters
uint64_t p, g, g_inv, beta, alpha;

//inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t c) {
//	return (uint64_t)(((__uint128_t)(a)*(__uint128_t)(b)) % c);
//}

//if no 128 bit registers are available use this function, but then we need n<32
inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t c) {
	return (uint64_t)(a*b) % c;
}

uint64_t modpow(uint64_t base, uint64_t exp, uint64_t modulus) {
	base %= modulus;
	uint64_t result = 1;
	while (exp > 0) {
		if (exp & 1) result = mulmod(result, base, modulus);
		base = mulmod(base, base, modulus);
		exp >>= 1;
	}

	return result;
}

uint64_t projection(uint64_t a)
{
	//get prg output based on a
	prg.seed(a);
	uint64_t val[] = { uint_dist(prg),uint_dist(prg),uint_dist(prg),uint_dist(prg),uint_dist(prg) };

	//create necessary masks
	uint64_t bitmask_prob = 1ULL << (max_digits);
	uint64_t bitmask_set = 1;
	uint64_t bitmask_prg = 1ULL << 63;
	int arr_count = 0;
	uint64_t res = 0;

	//create n bit vector, where each coordinate is one with current probability (weights[i-1]/n) based on the prg output
	for (int i = 0; i < n; ++i)
	{
		bool set = false;
		for (int j = 0; j < max_digits; ++j)
		{
			if (((val[arr_count] & bitmask_prg)>0) >((prob & bitmask_prob)>0))
				break;
			else if (((val[arr_count] & bitmask_prg)>0) < ((prob & bitmask_prob)>0))
			{
				set = true;
				break;
			}
			bitmask_prob >>= 1;
			bitmask_prg = (bitmask_prg == 1) ? 1ULL << 63 : bitmask_prg >> 1;
			if (bitmask_prg == (1ULL << 63)) arr_count++;
		}
		bitmask_prg = (bitmask_prg == 1) ? 1ULL << 63 : bitmask_prg >> 1;
		if (bitmask_prg == (1ULL << 63)) arr_count++;

		if (set)
			res ^= bitmask_set;
		bitmask_set <<= 1;
		bitmask_prob = 1ULL << max_digits;
	}

	//rejection sampling
	if ((res%p) == res)
		return res;
	return projection((a*flavour1 + flavour2) % p);
}
//compute inner product mod 2
bool scalar(uint64_t x, uint64_t y)
{
	uint64_t z = x&y;
	int wt = 0;
	while (z)
	{
		z &= (z - 1);
		wt += 1;
	}
	return (wt & 1) == 1;
}

//function subject to collision search
uint64_t function(uint64_t x)
{
	calls++;

	bool func = scalar(x, func_switch);
	x = (flavour1*x + flavour2) % p;
	x = projection(x);

	if (func)
		return modpow(g, x, p);
	return mulmod(modpow(g_inv, x, p), beta, p);
}

//Brent's collision finding algorithm
bool collision_finding(uint64_t sp)
{
	uint64_t slow = 1ULL << 63;
	uint64_t fast = function(sp);

	uint64_t  lam = 1, power = 1;
	//first phase: find colliding output
	while (slow != fast) {
		if (lam == power)
		{
			slow = fast;
			power <<= 1;
			lam = 0;
		}
		fast = function(fast);
		lam += 1;

	}

	//second phase: find colliding inputs
	slow = fast = sp;
	for (uint64_t i = 0; i<lam; ++i)
		fast = function(fast);
	if (sp == fast)
		return collision_finding((sp + 1) % p);

	uint64_t slow_tmp = 0, fast_tmp;
	while (slow != fast) {
		slow_tmp = slow;
		slow = function(slow);
		fast_tmp = fast;
		fast = function(fast);
	}

	coll[0] = slow_tmp;
	coll[1] = fast_tmp;

	return true;
}

//choose random number < p-1 with weight w
uint64_t rnd_num_with_wt(int w)
{
	int k = 0;
	uint64_t res = 0;
	while (k < w)
	{
		int r = uint_dist(mt) % n;
		if ((res&(1ULL << r)) == 0)
			res ^= (1ULL << r);
		k++;
		if (res >= (p - 1))
		{
			k = 0;
			res = 0;
		}
	}
	return res;
}
bool search_solution()
{
	do
	{
		//choose random flavours, startpoint and function switch
		func_switch = uint_dist(mt);
		if (n<64) func_switch %= bitmask_n;
		flavour1 = uint_dist(mt) % p;
		flavour2 = uint_dist(mt) % p;
		uint64_t sp = uint_dist(mt) % p;

		//start collision finding
		collision_finding(sp);

		//if collision in function one or function two, continue
		if (scalar(func_switch, coll[0]) == scalar(func_switch, coll[1]))
			continue;
		else {
			//we found the solution
			uint64_t x0 = (coll[0] * flavour1 + flavour2) % p;
			x0 = projection(x0);

			uint64_t x1 = (coll[1] * flavour1 + flavour2) % p;
			x1 = projection(x1);
			//final check
			if (((x0 + x1) % (p - 1)) == alpha)
				break;
		}

	} while (true);
	return true;
}

int main()
{
	n = 20;
	bitmask_n = 1ULL << n;

	//weight array at position i-1 states the weight needed for two numbers, such that they are expected to sum up to a number of weight i, in other words weights[i-1]=phi(i)*n
	//for n=20
	double weights[] = { 0.5070631622652736, 1.031825988341176, 1.5806462835255535, 2.161906885556572, 2.787443357625332, 3.475363040607318, 4.2564546365165885, 5.192112731632765, 6.447052936160102, 10.0 };
	//for n=32
	//double weights[] = { 0.5042188758684708, 1.018194848714795, 1.5440971419484837, 2.0844240683710975, 2.64212485564305, 3.2207703093044766, 3.8248016847891844, 4.459909372200531, 5.133637511363924, 5.856406460551055, 6.643371058093583, 7.518133665817576, 8.52121480956413, 9.73372557935704, 11.3730763392386, 16.0 };
	//for n=34
	//double weights[] = { 0.5039529833180001, 1.0169746564159423, 1.5409638893944821, 2.0780834351284567, 2.630850910283786, 3.2022628488893554, 3.7959699308007604, 4.416533728617788, 5.0698184042021595, 5.763616963079983, 6.508711014562079, 7.320797271375031, 8.224336922069577, 9.261335825258561, 10.515874496222104, 12.213233024112998, 17.0 };
	//for n=40
	//double weights[] = { 0.5033240704951129, 1.0141263245305472, 1.5337415083091333, 2.063651976682352, 2.605530085463872, 3.161292567051107, 3.733172079555371, 4.323813771113144, 4.936409100892026, 5.574886715250664, 6.244193710428533, 6.950726081214636, 7.703018038914196, 8.512909273033177, 9.397667082271969, 10.38422546326553, 11.518848729254078, 12.894105872320203, 14.757719370641526, 20.0 };
	//for n=48
	//double weights[] = { 0.5027419535981006, 1.0115379438228058, 1.5272922730721925, 2.0509874363013063, 2.583702603591918, 3.1266361025566463, 3.681133265564376, 4.248721530333432, 4.831155463956715, 5.43047557097868, 6.049086602534651, 6.689864058300797, 7.356302491977128, 8.052727665903724, 8.784609690826581, 9.559042667562355, 10.385513021412072, 11.277200498726364, 12.253342690529028, 13.343955815892773, 14.60058836903556, 16.126306692422492, 18.196499052216232, 24.0 };

	//p=previous_prime(1<<n), any other prime >2^(n-1) also possible
	if (n == 20)
		p = 1048573;
	if (n == 32)
		p = 4294967291;
	if (n == 34)
		p = 17179869143;
	if (n == 40)
		p = 1099511627689;	
	if (n == 48)
		p = 281474976710597;

	g = 2;
	g_inv = modpow(g, p - 2, p);

	//test performance

	int iterations = 10;
	uint64_t total = 0;
	for (int i = 2; i <= n / 2; ++i)
	{

		total = 0;
		int alpha_weight = i;

		//create max_digits-bit binary expansion of probability
		prob = ((1 << (max_digits + 1)) - 1)*1.0*weights[i - 1] / n;

		for (int j = 0; j < iterations; ++j)
		{
			//choose random dlog
			alpha = rnd_num_with_wt(alpha_weight);
			beta = modpow(g, alpha, p);

			//solve dlog
			search_solution();
			total += calls;
			calls = 0;
		}
		//output
		std::cout << "[" << i << "," << log2(total*1.0 / iterations) << "], ";
		fflush(stdout);
	}
	std::cout<<"\n";
}
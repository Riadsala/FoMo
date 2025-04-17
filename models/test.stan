transformed data {

vector[4] x = rep_vector(1, 4);
vector[4] z = rep_vector(1, 4);
vector[4] w;
x[2] = 2;
z[3] = 0;
z[4] = 0;

print("x: ", x);
print("z: ", z);

print("--------------------");

print("log x: ", log(x));
print("log z: ", log(z));

print("--------------------");

w = log(x)+log(z);
print("log x - log z: ", w);
print("log x .* z: ", log(x.*z));

print("--------------------");

print("log_sum_exp(w):", log_sum_exp(w));

print("--------------------");

w = w -log_sum_exp(w);
print("log normalised: ", w);
print("w:", exp(w));





}



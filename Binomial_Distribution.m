% Parameters
n = input('Enter the number of trials (n): ');  % Number of trials
p = input('Enter the probability of success (p): ');  % Probability of success
x = input('Enter the number of successes (x): ');  % Number of successes

% The probability mass function (pmf) of the binomial distribution
pmf = binopdf(x, n, p);

% The cumulative distribution function (cdf) of the binomial distribution
cdf = binocdf(x, n, p);

% The p-value
p_value = 1 - cdf;

% Display results
fprintf('Probability of getting exactly %d successes in %d trials: %.4f\n', x, n, pmf);
fprintf('Cumulative probability of getting up to %d successes in %d trials: %.4f\n', x, n, cdf);
fprintf('P-value of getting more than %d successes in %d trials: %.4f\n', x, n, p_value);

%%

% Define the parameters
x = 3; % Number of successes
y = 5;
% Number of trials
p = 0.1; % Probability of success on a single trial (e.g., for a fair coin)

% Calculate the binomial probability
binomial_probability = nchoosek(y, x) * (p^x) * ((1-p)^(y-x));

% Display the result
fprintf('The probability of getting exactly %d successes in %d trials is %.5f\n', x, y, binomial_probability);
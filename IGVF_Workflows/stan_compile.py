import cmdstanpy

models_to_compile = [
    '~/cmdstan/cs-guide-mixture.stan',
    '~/cmdstan/dc-guide-mixture.stan',
]

for model_path in models_to_compile:
    model = cmdstanpy.CmdStanModel(stan_file=model_path)
    model.compile()
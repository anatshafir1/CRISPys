import yaml

with open("configs_file.yml") as f:
    data = yaml.safe_load(f)


print(data)

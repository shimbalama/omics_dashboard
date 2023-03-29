def make_list_of_dicts(values: list[str]) -> list[dict[str, str]]:
    """Convert a list of strs into a list where those strings are values in dicts
    against keys label and value, for use in callbacks"""
    return [{"label": val, "value": val} for val in values]


''' Work out how many years of data the full MC count is equivalent to
'''
import json
from ConfigParser import ConfigParser, NoOptionError

def lt_equivs(rates, n_gens):
    lts = {}
    for name in rates.keys():
        rate  = rates[name]

        if n_gens[name] is None:
            lts[name] = -1
            continue
        else:
            count = n_gens[name]
            

        if rate != 0:
            lts[name] = count/rate
        else:
            lts[name] = -1

    return lts


def load_rates(event_config_file):
    ''' as dictionaries
    '''
    config = ConfigParser()
    config.read(event_config_file)
    
    sections = config.sections()
    rates = {}
    n_gens = {}
    for section in sections:
        if section == "summary":
            continue
        
        name = section
        rate = float(config.get(section, "rate"))
        try:
            n_gen = int(config.get(section, "n_generated"))
        except NoOptionError as e_:
            n_gen = None
        
        n_gens[section] = n_gen
        rates[section] = rate
    return rates, n_gens


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("event_config_file")
    parser.add_argument("output_file")
    parser.add_argument("-to_json", type=str, default=None, dest = "jsfile")

    args = parser.parse_args()
    
    rates, ngens = load_rates(args.event_config_file)
    ltes =  lt_equivs(rates, ngens)
    
    js_formatted = json.dumps(sorted(ltes.items(), key = lambda x : x[1], reverse = True), indent = 4, separators=(',', ': '))
    if args.jsfile is not None:
        with open(args.jsfile, "w") as f:
            f.write(js_formatted)

    with open(args.output_file, "w") as f:
        f.write("\n".join("{0: <40} | {1: <20}".format(x, y) for x, y in sorted(ltes.iteritems(), key = lambda x: x[1])))

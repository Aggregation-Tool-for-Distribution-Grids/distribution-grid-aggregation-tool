import pandapower as pp
import create_network_example

network = create_network_example.create_network()

# Export to Excel, JSON, Pickle
pp.to_excel(network, 'network.xlsx')
pp.to_json(network, 'network.json')
pp.to_pickle(network, 'network.p')

# Import Excel, JSON, Pickle
network = pp.from_excel('network.xlsx')
network = pp.from_json('network.json')
network = pp.from_pickle('network.p')



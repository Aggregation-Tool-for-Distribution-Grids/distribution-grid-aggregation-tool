import pandapower as pp

def create_network():

    # Empty network
    network = pp.create_empty_network()

    # Busses
    b_T = pp.create_bus(network, vn_kv=110.)
    b_150kV = pp.create_bus(network, vn_kv=110.)
    b_SUBS_A = pp.create_bus(network, vn_kv=20.)
    b_SUBS_A1 = pp.create_bus(network, vn_kv=20.)
    b_SUBS_A2 = pp.create_bus(network, vn_kv=20.)
    b_GEN_Y25 = pp.create_bus(network, vn_kv=0.4)
    b_GEN_Y54 = pp.create_bus(network, vn_kv=0.4)

    # External grid
    pp.create_ext_grid(network, bus=b_T, vm_pu=1.0, s_sc_max_mva=2182, s_sc_min_mva=598)

    # Lines
    pp.create_line(network, from_bus=b_T, to_bus=b_150kV, length_km=4.0, std_type="NAYY 4x50 SE")
    pp.create_line(network, from_bus=b_SUBS_A, to_bus=b_SUBS_A1, length_km=3.680, std_type="NAYY 4x50 SE")
    pp.create_line(network, from_bus=b_SUBS_A, to_bus=b_SUBS_A2, length_km=3.680, std_type="NAYY 4x50 SE")

    # Transformers
    pp.create_transformer(network, b_150kV, b_SUBS_A, std_type="25 MVA 110/20 kV", tap_pos=7)
    pp.create_transformer(network, b_SUBS_A1, b_GEN_Y25, std_type="0.63 MVA 20/0.4 kV")
    pp.create_transformer(network, b_SUBS_A2, b_GEN_Y54, std_type="0.63 MVA 20/0.4 kV")

    # Generators
    pp.create_sgen(network, b_GEN_Y25, -0.3, q_mvar=0, max_q_mvar=-0.1, min_q_mvar=0.1)
    pp.create_sgen(network, b_GEN_Y54, -0.3, q_mvar=0, max_q_mvar=-0.05, min_q_mvar=0.05)

    # Switches
    pp.create_switch(network, b_SUBS_A1, b_SUBS_A2, "b", closed=False)

    return network
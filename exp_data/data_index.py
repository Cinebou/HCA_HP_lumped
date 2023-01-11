# data index of the file, after 2022/12/21
def index_data_ver3():
    param_d = {
            "time"                  : "日付/時間",
            # flow rate
            "MF1-CO2 flow_"          : "CH1",  
            "MF2-Water1_"            : "CH2",  
            "MF3-Water2"            : "CH3",

            # pressure  
            "PS1-compressor inlet"  : "CH6",
            "PS2-compressor outlet" : "CH7",
            "Inlet"                 : "CH8",
            "Outlet"                : "CH9",
            "PS5-tank2 expV"        : "CH10",
            "PS6-tank2 comp"        : "CH11",
            "PS7-WHX1 comp"         : "CH12",
            "PS8-WHX1 expV"         : "CH13",
            "PS9-WHX2 expV"         : "CH14",
            "PS10-WHX2 comp"        : "CH15",

            # temp
            #"T-room"                : "CH20",
            "T1-compressor inlet"   : "CH21",
            "T2-compressor outlet"  : "CH22",
            "T3-blank"              : "CH23",
            "T4-blank"              : "CH24",
            "T5-tank2 expV"         : "CH25",
            "T6-tank2 comp"         : "CH26",
            "T7-WHX1 comp"          : "CH27",
            "T8-WHX1 expV"          : "CH28",
            "T9-WHX2 expV"          : "CH29",
            "T10-WHX2 comp"         : "CH30",
            "T11-tank1 w inlet"     : "CH41",
            "T12-tank1 w outlet"    : "CH42",
            "T13-tank2 water inlet" : "CH43",
            "T14-tank2 water outlet": "CH44",
            "T21-tank1 mof"         : "CH45",
            "T22-tank1 mof"         : "CH46",
            "T23-tank1 mof"         : "CH47",
            "T24-tank1 gas low"     : "CH48",

            # sheath thermocouples
            "T25-co2 inlet"         : "CH49",
            "T26-co2 outlet"        : "CH50",

            # tank2 thermocouples
            "T27-tank2 mof"         : "CH51",
            "T28-tank2 mof"         : "CH52",
            "T29-tank1 gas high"    : "CH53",
            "T30-tank2 gas"         : "CH54",

            # vessel temp
            "T51-tank1 vessel CO2in" : "CH61",
            "T52-tank1 vessel bottom": "CH62",
            "T53-tank1 vessel bottom": "CH63",
            "T54-tank1 vessel bottom": "CH64",
            "T55-tank1 vessel bottom": "CH65",
            "T56-tank1 vessel middle": "CH66",
            "T57-tank1 vessel middle": "CH67",
            "T58-tank1 vessel middle": "CH68",
            "T59-tank1 vessel middle": "CH69",
            "T60-tank1 vessel upper" : "CH70",
            "T61-tank1 vessel upper" : "CH71",    
            "T62-tank1 vessel upper" : "CH72",        
            "T63-tank1 vessel upper" : "CH73",
            "T64-tank1 vessel top"   : "CH74",

            # room temp
            "T-room"               : "CH75",

            # circuit temp
            "T71 circuit NV4"        : "CH81",
            "T72 circuit before BV1" : "CH82",
            "T73 circuit after BV1"  : "CH83",
            "T74 circuit before TV1" : "CH84",
            "T75 circuit before expV": "CH85",
            "T76 circuit after expV" : "CH86",
            "T77 circuit after TV4"  : "CH87",
            "T78 circuit before BV3" : "CH88",
            "T79 circuit before MS1" : "CH89",
            "T80 circuit after MS1"  : "CH90"

    }
    return param_d
import mesa
from cancerInmunoediting.model import CancerInmunoediting
# dictionary of user settable parameters - these map to the model __init__ parameters
model_params = {
    "meanIS": mesa.visualization.Slider(
        "meanIS", 0.9, 0.01, 1, 0.01, description="Mean of the Inmune System distribution"
    ),
    "stdIS": mesa.visualization.Slider(
        "stdIS", 0.05, 0.01, 1, 0.01, description="Standard Deviation of the Inmune System distribution",
    ),
    "meanCancer": mesa.visualization.Slider(
        "meanCancer", 0.9, 0.01, 1, 0.01, description="Mean of the Cancer distribution",
    ),
    "stdCancer": mesa.visualization.Slider(
        "stdCancer", 0.05, 0.01, 1, 0.01, description="Standard deviation of the Cancer distribution",
    ),
    "ml": mesa.visualization.Slider(
        "ml", 100, 0.001, 1000, 0.01, description="Standard deviation of the Cancer distribution",
    ),
}
# map data to chart in the ChartModule
chart_element1 = mesa.visualization.ChartModule(
    [
        {"Label": "CancerCells", "Color": "#2596be"},
        {"Label": "M1",     "Color": "#9925be"},
        {"Label": "M2",     "Color": "#be4d25"},
        {"Label": "N1",     "Color": "#49be25"},
        {"Label": "N2",     "Color": "#bea925"},
        {"Label": "NK",     "Color": "#041014"},
        {"Label": "T",     "Color": "#873e23"},
        {"Label": "Th",     "Color": "#F90BB6"},
        {"Label": "Treg",     "Color": "#F9200B"},
    ]
)
chart_element2 = mesa.visualization.ChartModule(
    [
        {"Label": "ProCancer",  "Color": "#2596be"},
        {"Label": "AntiCancer", "Color": "#041014"}
    ]
)

chart_element3 = mesa.visualization.ChartModule(
    [
        {"Label": "HAntiCancer",  "Color": "#2596be"},
        {"Label": "HProCancer", "Color": "#041014"},
        {"Label": "HTME", "Color": "#be4d25"},
    ]
)
chart_element4 = mesa.visualization.ChartModule(
    [
        {"Label": "TumorGrowthRate",  "Color": "#2596be"},
    ]
)
# create instance of Mesa ModularServer
server = mesa.visualization.ModularServer(
    CancerInmunoediting,
    [chart_element1, chart_element2, chart_element3,chart_element4],
    "Cancer Inmunoediting Model",
    model_params=model_params,
)

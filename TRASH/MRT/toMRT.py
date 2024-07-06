from astropy.table import Table
import cdspyreadme

tablemaker = cdspyreadme.CDSTablesMaker()
tablemaker.title = "GAMAV"
tablemaker.author = 'O.Ryzhov'

csv = Table.read('MRT_my.csv')
#
colcataid = csv["CATAID"]
colcataid.name = "GAMA_ID"
colcataid.description="Unique GAMA object ID"
colcataid.unit='-'

colra = csv["RA"]
colra.name = "RA"
colra.description="J2000"
colra.unit='deg'

coldec = csv["DEC"]
coldec.name = "DEC"
coldec.description="J2000"
coldec.unit='deg'

colz = csv["Z"]
colz.name = "Z" #to check
colz.description="Redshift with heliocentric correction removed"
colz.unit='-'

colbms = csv["bMS/MS"]
colbms.name = "bMS_MS"
colbms.description="Galaxy position with the respect to the Main Sequence (eq. 28, Speagle et al. 2014). Below-MS: 0, MS: 1"
colbms.unit='-'

colbpt = csv["BPT"]
colbpt.name = "BPT"
colbpt.description="Spectral classification via BPT diagram (classes described in Ryzhov et al. 2024)"
colbpt.unit='-'

colwhan = csv["WHAN"]
colwhan.name = "WHAN"
colwhan.description="Spectral classification via WHAN diagram (classes described in Ryzhov et al. 2024)"
colwhan.unit='-'

colon16 = csv["outflow_agn_on_percentile16"]
colon16.name = "outflow_agn_on_percentile16"
colon16.description="16th percentile of outflow of the molecular gas concerning AGN luminosity"
colon16.unit='dex(mass/year)'

colon50 = csv["outflow_agn_on_percentile50"]
colon50.name = "outflow_agn_on_percentile50"
colon50.description="50th percentile of outflow of the molecular gas concerning AGN luminosity"
colon50.unit='dex(mass/year)'

colon84 = csv["outflow_agn_on_percentile84"]
colon84.name = "outflow_agn_on_percentile84"
colon84.description="84th percentile of outflow of the molecular gas concerning AGN luminosity"
colon84.unit='dex(mass/year)'

coloff16 = csv["outflow_agn_off_percentile16"]
coloff16.name = "outflow_agn_off_percentile16"
coloff16.description="16th percentile of outflow of the molecular gas without concerning AGN luminosity"
coloff16.unit='dex(mass/year)'

coloff50 = csv["outflow_agn_off_percentile50"]
coloff50.name = "outflow_agn_off_percentile50"
coloff50.description="50th percentile of outflow of the molecular gas without concerning AGN luminosity"
coloff50.unit='dex(mass/year)'

coloff84 = csv["outflow_agn_off_percentile84"]
coloff84.name = "outflow_agn_off_percentile84"
coloff84.description="84th percentile of outflow of the molecular gas without concerning AGN luminosity"
coloff84.unit='dex(mass/year)'

table = tablemaker.addTable('MRT_my.csv', description='My catalogue')
tablemaker.toMRT()

tablemaker.makeReadMe()

with open("ReadMe", "w") as fd:
    tablemaker.makeReadMe(out=fd)
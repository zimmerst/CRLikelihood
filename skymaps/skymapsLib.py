# @file skymapsLib.py
# @brief scons package dependencies for skymaps
#
#$Header: /nfs/slac/g/glast/ground/cvs/skymaps/skymapsLib.py,v 1.8 2010/12/09 00:01:46 kerrm Exp $
def generate(env, **kw):
	if not kw.get('depsOnly',0):
		env.Tool('addLibrary', library=['skymaps'])
        depends = 'facilities tip st_facilities healpix embed_python timeSystem'.split()
        for pack in depends: env.Tool(pack+'Lib')

	# why these?
	#env.Tool('addLibrary', library=env['clhepLibs'])
	#env.Tool('addLibrary', library=env['cfitsioLibs'])

def exists(env):
	return 1

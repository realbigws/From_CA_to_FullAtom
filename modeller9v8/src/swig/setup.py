from distutils.core import setup, Extension
import commands

# Change this to your Modeller executable type:
exetype = "i386-intel8"

def pkgconfig(*packages, **kw):
    """Utility function to parse pkg-config output"""
    flag_map = {'-I': 'include_dirs', '-L': 'library_dirs', '-l': 'libraries'}
    for token in commands.getoutput("pkg-config --libs --cflags %s" \
                                    % ' '.join(packages)).split():
        kw.setdefault(flag_map.get(token[:2]), []).append(token[2:])
    return kw

# Get paths for glib 2.0:
glib = pkgconfig("glib-2.0")

mod = Extension('_modeller', sources=['modeller_wrap.c'],
                include_dirs=['../include', '../include/%s' % exetype] + \
                             glib['include_dirs'],
                libraries=['modeller'] + glib['libraries'],
                library_dirs=['../../lib/%s' % exetype] \
                             + glib.get('library_dirs', []))

setup(name='Modeller',
      description='Protein structure modeling by satisfaction of ' \
                  + 'spatial restraints',
      author='Andrej Sali', url='http://www.salilab.org/modeller/',
      ext_modules=[mod])

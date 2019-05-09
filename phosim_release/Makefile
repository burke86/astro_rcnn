all:
	 rm -f source/raytrace/*.o; rm -f source/trim/*.o; rm -f source/e2adc/*.o; rm -f source/instrument/*.o; rm -f source/atmosphere/*.o ; rm -f source/ancillary/*.o; rm -f source/catalogtools/*.o
	 csh -c "grep -q EXT bin/setup && echo '#define EXT' >! source/raytrace/ext.h || echo '#undef EXT' >! source/raytrace/ext.h"
	 csh -c "source bin/setup ; cd source/raytrace ; make; mv raytrace ../../bin/."
	 csh -c "source bin/setup ; cd source/trim ; make; mv trim ../../bin/."
	 csh -c "source bin/setup ; cd source/e2adc ; make; mv e2adc ../../bin/."
	 csh -c "source bin/setup ; cd source/instrument ; make; mv instrument ../../bin/."
	 csh -c "source bin/setup ; cd source/atmosphere ; make; mv atmosphere ../../bin/."
	 csh -c "source bin/setup ; cd source/ancillary ; make "
	 csh -c "source bin/setup ; cd source/catalogtools ; make;  mv phosimcatgen  ../../bin/."
	 csh -c "source bin/setup ; cd validation ; make "
	 csh -c "rm -f bin/version"
	 csh -c "if (-d .git) git log -1 --format=commit' '%H >! bin/version"
	 csh -c "if (-d .git) git describe --tags >> bin/version"
	 csh -c "if (! -d .git) echo 'commit none' >> bin/version"
	 csh -c "if (! -d .git) pwd >> bin/version"
	 @echo ' '
	 @echo 'The Photon Simulator has been compiled successfully!'
	 @echo ' '
	 @echo 'Type "phosim_gui" or "phosim" to begin!'
	 @echo ' '

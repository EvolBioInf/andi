#
# Regular cron jobs for the andi package
#
0 4	* * *	root	[ -x /usr/bin/andi_maintenance ] && /usr/bin/andi_maintenance

/** Prototype for callback functions accepted by set_log_function() */
typedef int (*cb_log)(void *logdata, const char *text);

/** Write text directly to the log */
void mod_log_write(const char *str);

/** Get the current log level */
int mod_log_get(int level);

/** Set the current log level */
void mod_log_set(int level, int value);

/** Install a callback function to handle all log output */
void mod_log_function_set(cb_log logfunc, cb_free freefunc, void *logdata);

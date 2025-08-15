import pyodbc
import traceback
import os
import pandas as pd
from common_functions import *


class Database:
    def __init__(self, db_cfg):
        
        # load configuration values
        self.config = db_cfg
        self.s_conn = self.prepare_conn_string()
        self.conn = None

    def prepare_conn_string(self):
        # get connection string template
        # print (self.config)

        conn_str = get_config_value('mdb_conn_str', self.config).strip()
        
        # env_db_server = get_config_value('database/connection/env_db_server', self.config).strip()
        # print('env_db_server name = {}'.format(env_db_server))
        # print('env_db_server value = {}'.format(os.environ.get(env_db_server)))
        
        # get values for connection string components
        driver = get_and_validate_env_var(get_config_value('env_db_driver', self.config).strip())
        server = get_and_validate_env_var(get_config_value('env_db_server', self.config).strip())
        dbname = get_and_validate_env_var(get_config_value('env_db_name', self.config).strip())
        user_name = get_and_validate_env_var(get_config_value('env_db_user_name', self.config).strip())
        user_pwd = get_and_validate_env_var(get_config_value('env_db_user_pwd', self.config).strip())
        
        db_driver_plh = get_config_value('db_driver_plh', self.config).strip()
        server_plh = get_config_value('server_plh', self.config).strip()  # '|server|'
        dbname_plh = get_config_value('dbname_plh', self.config).strip()  # '|db_name|'
        user_name_plh = get_config_value('user_name_plh', self.config).strip()  # '|db_user_name|'
        user_pwd_plh = get_config_value('user_pwd_plh', self.config).strip()  # '|db_user_pwd|'
        
        # update connection string template with retrieved values
        conn_str = conn_str.replace(db_driver_plh, driver)
        conn_str = conn_str.replace(server_plh, server)
        conn_str = conn_str.replace(dbname_plh, dbname)
        conn_str = conn_str.replace(user_name_plh, user_name)
        conn_str = conn_str.replace(user_pwd_plh, user_pwd)
        # print('conn_str = {}'.format(conn_str))  # for testing only
        return conn_str

    def open_connection(self):
        print('Attempting to open connection to Metadata DB.')
        try:
            self.conn = pyodbc.connect(self.s_conn, autocommit=True)
            print('Successfully established the database connection.')
        except Exception as ex:
            # report unexpected error during openning DB connection
            _str = 'Unexpected Error "{}" occurred during an attempt to open connecton to database ({})\n{} ' \
                .format(ex, self.s_conn, traceback.format_exc())
            print(_str)
            raise


class DBAccess(Database):
    def __init__(self, db_cfg):
        Database.__init__(self, db_cfg)

    def execute_query_with_output(self, str_sql):  # sample_id, row_json, dict_json, filepath):
        if not self.conn:
            # open connection if needed
            self.open_connection()

        if not self.conn:
            # connection failed to open
            print('Database connection cannot be established (see earlier errors for details).')
            return None

        print('Attempting to execute the following SQL call: {}'.format(str_sql.rstrip()))

        try:
            # get reference to connection object
            cursor = self.conn.cursor()
            
            # execute stored procedure and get output to pandas dataframe
            
            df = pd.read_sql(str_sql, self.conn)
            return df

        except Exception as ex:
            # report an error if DB call has failed.
            _str = 'Error "{}" occurred during execution of the following SQL script "{}". ' \
                   'Here is the traceback: \n{} '.format(
                ex, str_sql, traceback.format_exc())
            print(_str)
            # raise the error further to break execution of the snakemake
            raise (ex)
            



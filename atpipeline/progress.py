'''
Created on Dec 8, 2012

@author: uemit.seren
'''

import pika
import simplejson


class MessengerBase:

    def __init__(self):
        self.progress = 0.0
        self.step = 0.01
        self.task = ''
        
    def update_progress_bar(self, progress=None, step=None, task_status=''):
        self.update_status(progress, step, task_status)
    
    def update_status(self, progress=None, step=None, task_status=''):
        if progress is not None:
            self.progress = progress
        elif step is not None:
            self.progress = self.progress + step
        else:
            self.progress = self.progress + self.step
        if self.progress > 1.0:
            self.progress = 1.0
        self.task_status = task_status
        self.send_message()

    
    def set_step(self, step):
        self.step = step

    def close_file(self):
        return
        

class StdoutMessenger(MessengerBase):
    
    def send_message(self):
        print "Progress: %s : Status %s" % (self.progress,self.task_status)

        

class ProgressMessenger(MessengerBase):


    def __init__(self, host, port, username, password):
        MessengerBase.__init__(self)
        credentials = pika.PlainCredentials(username, password)
        parameters = pika.ConnectionParameters(host, port, '/', credentials)
        self.host = host
        self.connection = pika.BlockingConnection(parameters)
        self.channel = self.connection.channel()
        self.channel.queue_declare(queue='gwaspipeline')
    
    def send_message(self):
        body = {'progress':self.progress,'task_status':self.task_status}
        self.channel.basic_publish(exchange='',routing_key='gwaspipeline',body=simplejson.dumps(body))

    def close_file(self):
        self.connection.close()
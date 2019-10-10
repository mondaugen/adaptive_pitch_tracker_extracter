gdb.execute('source src/test/common.gdb')
gdb.execute('file %s' % ("src/test/bin/stretch_f32_file",))
ex=gdb.execute

class Breakpoint_dft_forward(gdb.Breakpoint):
    def __init__(self,exp,window_length):
        gdb.Breakpoint.__init__(self,exp)
        self.window_length=window_length
        print(self.window_length)
        self.addr_keys=[]

    def stop(self):
        # arguments are expression-giving-pointer length output-file
        addr=int(gdb.parse_and_eval('b'))
        if addr in self.addr_keys:
            cmd='a_rz64'
        else:
            cmd='d_rz64'
            self.addr_keys.append(addr)
        ex('%s a %d /tmp/%u.z64' % (cmd,self.window_length,int(addr)))
        return False

class Breakpoint_pvs_process(gdb.Breakpoint):
    def stop(self):
        # get window length
        window_length=gdb.parse_and_eval('pvs->config.window_length')
        # print addresses of the stuff we will dump
        print('z_input0 %u' % (int(gdb.parse_and_eval('pvs->z_input0'))))
        print('z_inputH %u' % (int(gdb.parse_and_eval('pvs->z_inputH'))))
        print('z_outputH %u' % (int(gdb.parse_and_eval('pvs->z_outputH'))))
        # have to add offset address because the arguments aren't valid yet
        self.dft_forward_break=Breakpoint_dft_forward(
            '*pvs->func_table->math.dft_forward + 24',
            int(window_length))
        self.enabled=False
        return False
        

b1=Breakpoint_pvs_process('pvs_process')

gdb.execute('run')

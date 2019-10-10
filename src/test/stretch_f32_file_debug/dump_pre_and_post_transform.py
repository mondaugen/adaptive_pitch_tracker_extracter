gdb.execute('source src/test/common.gdb')
gdb.execute('file %s' % ("src/test/bin/stretch_f32_file",))
ex=gdb.execute

class post_DFT_break(gdb.FinishBreakpoint):
    def __init__(self,super_break,b_addr,a_addr):
        gdb.FinishBreakpoint.__init__(self,internal=True)
        self.super_break=super_break
        self.b_addr=b_addr
        self.a_addr=a_addr
    def stop(self):
        # this dumps the contents of the frequency domain buffer (e.g.,
        # z_input0) after the transform     
        if self.b_addr in self.super_break.addr_keys:
            cmd='a_rz64'
            r_cmd='a_f32'
        else:
            cmd='d_rz64'
            r_cmd='d_f32'
            self.super_break.addr_keys.append(self.b_addr)
        scmd=('%s %u %d /tmp/%u.z64' % 
            (cmd,self.b_addr,self.super_break.window_length,self.b_addr))
        #print(scmd)
        ex(scmd)
        # this dumps the contents from the time domain buffer (i.e., r_workspace)
        # after the transform
        # we identify it with the address of the output buffer so we can see
        # where the corresponding output buffer came from 
        scmd=('%s %u %d /tmp/%u.f32' % 
            (r_cmd,self.a_addr,self.super_break.window_length,self.b_addr))
        #print(scmd)
        ex(scmd)
        return False

class Breakpoint_dft_forward(gdb.Breakpoint):
    def __init__(self,exp,window_length):
        gdb.Breakpoint.__init__(self,exp)
        self.window_length=window_length
        self.addr_keys=[]
        self.cur_post_break=None

    def stop(self):
        # the address of the b function parameter, which is the output where the
        # frequency domain is put for dft_forward
        b_addr=int(gdb.parse_and_eval('b'))
        # the address of the a function parameter, which is the output where the
        # time domain is put for dft_forward
        a_addr=int(gdb.parse_and_eval('a'))
        self.cur_post_break=post_DFT_break(self,b_addr,a_addr)
        return False

class Breakpoint_pvs_process(gdb.Breakpoint):
    def stop(self):
        # get window length
        window_length=gdb.parse_and_eval('pvs->config.window_length')
        with open('/tmp/dump_pre_and_post_transform.vars','w') as f:
            f.write("W=%d\n" % (window_length,))
            # print addresses of the stuff we will dump
            f.write('z_input0=%u\n' % (int(gdb.parse_and_eval('pvs->z_input0'))))
            f.write('z_inputH=%u\n' % (int(gdb.parse_and_eval('pvs->z_inputH'))))
            f.write('z_outputH=%u\n' % (int(gdb.parse_and_eval('pvs->z_outputH'))))
        # have to add offset address because the arguments aren't valid yet
        self.dft_forward_break=Breakpoint_dft_forward(
            '*pvs->func_table->math.dft_forward + 24',
            int(window_length))
        self.enabled=False
        return False
        

b1=Breakpoint_pvs_process('pvs_process')

gdb.execute('run')

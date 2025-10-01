# import:
import qiskit
from qiskit import QuantumCircuit
from qiskit.circuit import QuantumRegister, ClassicalRegister
from qiskit_aer import AerSimulator
import math

def w(n):
    cnt = 0
    while(n > 0):
        cnt += (n & 1)
        n = n >> 1
    return cnt

def r(n, t):
    # 상수 m의 값을 n과 t로부터 계산하는 함수.
    return int(math.floor(n / (2 ** t)))

def GetBitLen(num1):
    result = 0
    if(num1 < 2):
        result = 1
    else:
        result = int(math.floor(math.log2(num1))) + 1
    return result

def Adder_Draper(op1, op2, carry, ancilla, targetCircuit, zerorize_carry=True, is_op2_qubit=True):
    # Drapper's adder
    '''
    Adder_Draper(op1, op2, carry, ancilla, targetCircuit, zerorize_carry, is_op2_qubit)
    Do:
    - op1 = op1 + op2
    on targetCircuit
    Operands type:
    - op1: QuantumRegister
    - op2: QuantumRegister or integer
    - carry: QuantumRegister
    - ancilla: QuantumRegister
    - targetCircuit: QuantumCircuit
    - zerorize_carry: boolean
    - is_op2_qubit: boolean
    Description:
    - this fuction is Draper's in-place adder.
    - the paper can be accessed with:
    https://arxiv.org/abs/quant-ph/0406142
    or can be found with:
    T. G. Draper, S. A. Kutin, E. M. Rains, and K. M. Svore, “A logarithmic-depth quantum carry-lookahead adder,” 2004.
    - if zerorize_carry is True, then carry qubit will be zerorized.
    - if zerorize_carry is False, then carry qubit will not be zerorized and store carry of the addition.
    - if is_op2_qubit is True, then op2 should be qubit.
    - if is_op2_qubit is False, then op2 should be classical integer value.
    '''
    #
    # set parameters:
    z = carry # short name of carry qubit
    x = ancilla # short name of ancilla qubit
    n = 0 # n of "n-bit addition"
    if is_op2_qubit:
        n = len(op2)
    else:
        n = GetBitLen(op2)
    n = max(n, len(op1))

    # preprocess if op2 is classic value:
    op2_list = None # listfied each bit of op2
    if not is_op2_qubit:
        op2_list = list(map(int, bin(op2)[2:].zfill(len(op1))))[::-1]
    #

    # 1:
    loopAmount = 0
    if zerorize_carry:
        loopAmount = n-1
    else:
        loopAmount = n
    for i in range(0, loopAmount, 1):
        if is_op2_qubit:
            targetCircuit.ccx(op1[i], op2[i], z[i]) # z[i] = op1[i] * op2[i] = carry[i+1]
        elif op2_list[i] == 1:
            targetCircuit.cx(op1[i], z[i]) # z[i] = op1[i] * op2[i] = carry[i+1]
    #
    # 2:
    loopAmount = n
    for i in range(0, loopAmount, 1):
        if is_op2_qubit:
            targetCircuit.cx(op2[i], op1[i]) # op1[i] = op1[i] ^ op2[i] = sum[i]
        elif op2_list[i] == 1:
            targetCircuit.x(op1[i]) # op1[i] = op1[i] ^ op2[i] = sum[i]
    #
    # 3:
    # P-round:
    idx_x = 0
    idx_x_pre_t = 0
    loopAmount = int(math.floor(math.log2(n)))
    for t in range(1, loopAmount, 1):
        idx_x_pre_t = idx_x - (r(n, (t - 1)) - 1)
        for m in range(1, r(n, t), 1):
            if(t == 1):
                targetCircuit.ccx(op1[2 * m], op1[(2 * m) + 1], x[idx_x])
                #print("({}, {}) | x[{}] ^= op1[{}] * op1[{}]".format(t, m, idx_x, (2 * m), (2 * m + 1)))
            else:
                targetCircuit.ccx(x[idx_x_pre_t - 1 + (2 * m)], x[idx_x_pre_t + (2 * m)], x[idx_x])
                #print("({}, {}) | x[{}] ^= x[{}] * x[{}]".format(t, m, idx_x, (idx_x_pre_t - 1 + (2 * m)), (idx_x_pre_t + (2 * m))))
            idx_x += 1
    #
    # G-round:
    idx_x = 0
    idx_x_pre_t = 0
    loopAmount = int(math.floor(math.log2(n))) + 1
    for t in range(1, loopAmount, 1):
        idx_x_pre_t = idx_x - (r(n, (t - 1)) - 1)
        for m in range(0, r(n, t), 1):
            if zerorize_carry:
                if((t == 1) and (((2 * m) + 1) != (n - 1))):
                    targetCircuit.ccx(z[(2 * m)], op1[(2 * m) + 1], z[(2 * m) + 1])
                    #print("({}, {}) | z[{}] ^= z[{}] * op1[{}]".format(t, m, (2 * m + 1), (2 * m), (2 * m + 1)))                  
                elif((t != 1) and ((((2 ** t) * m) + (2 ** t) - 1) != (n - 1))):
                    targetCircuit.ccx(z[((2 ** t) * m) + (2 ** (t - 1)) - 1], x[(idx_x_pre_t + (2 * m))], z[((2 ** t) * m) + (2 ** t) - 1])
                    #print("({}, {}) | z[{}] ^= z[{}] * x[{}]".format(t, m, (((2 ** t) * m) + (2 ** t) - 1), (((2 ** t) * m) + (2 ** (t - 1)) - 1), (idx_x_pre_t + (2 * m))))
            else:
                if(t == 1):
                    targetCircuit.ccx(z[(2 * m)], op1[(2 * m) + 1], z[(2 * m) + 1])
                    #print("({}, {}) | z[{}] ^= z[{}] * op1[{}]".format(t, m, (2 * m + 1), (2 * m), (2 * m + 1)))
                else:
                    targetCircuit.ccx(z[((2 ** t) * m) + (2 ** (t - 1)) - 1], x[(idx_x_pre_t + (2 * m))], z[((2 ** t) * m) + (2 ** t) - 1])
                    #print("({}, {}) | z[{}] ^= z[{}] * x[{}]".format(t, m, (((2 ** t) * m) + (2 ** t) - 1), (((2 ** t) * m) + (2 ** (t - 1)) - 1), (idx_x_pre_t + (2 * m))))
            idx_x += 1
        idx_x -= 1
    #
    # C-round:
    idx_x = 0
    loopAmount = int((math.floor(math.log2(n / 3))) + 1)
    for t in range(1, loopAmount, 1):
        idx_x += (r(n, t) - 1)
    idx_x_pre_t = 0
    for t in range(loopAmount, 0, -1):
        idx_x -= (r(n, (t - 1)) - 1)
        idx_x_pre_t = idx_x
        for m in range(1, r((n - (2 ** (t - 1))), t) + 1, 1):
            if zerorize_carry:
                if((t == 1) and ((2 * m) != (n - 1))):
                    targetCircuit.ccx(z[2 * m - 1], op1[2 * m], z[(2 * m)])
                    #print("({}, {}) | z[{}] ^= z[{}] * op1[{}]".format(t, m, (2 * m), (2 * m - 1), (2 * m)))
                elif((t != 1) and ((((2 ** t) * m) + (2 ** (t - 1)) - 1)  != (n - 1))):
                    targetCircuit.ccx(z[(2 ** t)  * m - 1], x[(idx_x_pre_t - 1 + (2 * m))], z[((2 ** t) * m) + (2 ** (t - 1)) - 1])
                    #print("({}, {}) | z[{}] ^= z[{}] * x[{}]".format(t, m, (((2 ** t) * m) + (2 ** (t - 1)) - 1), ((2 ** t)  * m - 1), (idx_x_pre_t - 1 + (2 * m))))
            else:
                if(t == 1):
                    targetCircuit.ccx(z[2 * m - 1], op1[2 * m], z[(2 * m)])
                    #print("({}, {}) | z[{}] ^= z[{}] * op1[{}]".format(t, m, (2 * m), (2 * m - 1), (2 * m)))
                else:
                    targetCircuit.ccx(z[(2 ** t)  * m - 1], x[(idx_x_pre_t - 1 + (2 * m))], z[((2 ** t) * m) + (2 ** (t - 1)) - 1])
                    #print("({}, {}) | z[{}] ^= z[{}] * x[{}]".format(t, m, (((2 ** t) * m) + (2 ** (t - 1)) - 1), ((2 ** t)  * m - 1), (idx_x_pre_t - 1 + (2 * m))))
    #
    # inverse P-round:
    idx_x = 0
    loopAmount = int(math.floor(math.log2(n))) - 1
    for t in range(1, loopAmount, 1):
        idx_x += (r(n, t) - 1)
    idx_x_pre_t = 0
    for t in range(loopAmount, 0, -1):
        idx_x_pre_t = idx_x - (r(n, (t - 1)) - 1)
        for m in range(1, r(n, t), 1):
            if(t == 1):
                targetCircuit.ccx(op1[2 * m], op1[(2 * m) + 1], x[idx_x])
                #print("({}, {}) | x[{}] ^= op1[{}] * op1[{}]".format(t, m, idx_x, (2 * m), (2 * m + 1)))
            else:
                targetCircuit.ccx(x[idx_x_pre_t - 1 + (2 * m)], x[idx_x_pre_t + (2 * m)], x[idx_x])
                #print("({}, {}) | x[{}] ^= x[{}] * x[{}]".format(t, m, idx_x, (idx_x_pre_t - 1 + (2 * m)), (idx_x_pre_t + (2 * m))))
            idx_x += 1
        idx_x = idx_x_pre_t
    #
    # 4:
    loopAmount = n
    for i in range(1, loopAmount, 1):
        targetCircuit.cx(z[i - 1], op1[i])
    # 5:
    loopAmount = n - 1
    for i in range(0, loopAmount, 1):
        targetCircuit.x(op1[i])
    #
    # 6:
    loopAmount = n - 1
    for i in range(1, loopAmount, 1):
        if is_op2_qubit:
            targetCircuit.cx(op2[i], op1[i])
        elif op2_list[i] == 1:
            targetCircuit.x(op1[i])
    #
    # 7:
    # reverse: inverse P round:
    idx_x = 0
    idx_x_pre_t = 0
    loopAmount = int(math.floor(math.log2(n)))
    for t in range(1, loopAmount, 1):
        idx_x_pre_t = idx_x - (r(n, (t - 1)) - 1)
        for m in range(1, r(n, t), 1):
            if((t == 1) and ((2 * m + 1) != (n - 1))):
                targetCircuit.ccx(op1[2 * m], op1[(2 * m) + 1], x[idx_x])
                #print("({}, {}) | x[{}] ^= op1[{}] * op1[{}]".format(t, m, idx_x, (2 * m), (2 * m + 1)))
            elif((t != 1)):
                targetCircuit.ccx(x[idx_x_pre_t - 1 + (2 * m)], x[idx_x_pre_t + (2 * m)], x[idx_x])
                #print("({}, {}) | x[{}] ^= x[{}] * x[{}]".format(t, m, idx_x, (idx_x_pre_t - 1 + (2 * m)), (idx_x_pre_t + (2 * m))))
            idx_x += 1
    #
    # reverse: C-round:
    idx_x = 0
    idx_x_pre_t = -(r(n, 0) - 1)
    loopAmount = int(math.floor(math.log2(n / 3))) + 2
    for t in range(1, loopAmount, 1):
        for m in range(r((n - (2 ** (t - 1))), t), 0, -1):
            if((t == 1) and ((2 * m) != (n - 1))):
                targetCircuit.ccx(z[2 * m - 1], op1[2 * m], z[(2 * m)])
                #print("({}, {}) | z[{}] ^= z[{}] * op1[{}]".format(t, m, (2 * m), (2 * m - 1), (2 * m)))
            elif((t != 1) and ((((2 ** t) * m) + (2 ** (t - 1)) - 1)  != (n - 1))):
                targetCircuit.ccx(z[(2 ** t)  * m - 1], x[(idx_x_pre_t - 1 + (2 * m))], z[((2 ** t) * m) + (2 ** (t - 1)) - 1])
                #print("({}, {}) | z[{}] ^= z[{}] * x[{}]".format(t, m, (((2 ** t) * m) + (2 ** (t - 1)) - 1), ((2 ** t)  * m - 1), (idx_x_pre_t - 1 + (2 * m))))
            '''
            if((t == 1) and ((2 * m) == (n - 1))):
                print("({}, {}) | z[{}] ^= z[{}] * op1[{}] | STRIPPED".format(t, m, (2 * m), (2 * m - 1), (2 * m)))
            elif((t != 1) and ((((2 ** t) * m) + (2 ** (t - 1)) - 1)  == (n - 1))):
                print("({}, {}) | z[{}] ^= z[{}] * x[{}] | STRIPPED".format(t, m, (((2 ** t) * m) + (2 ** (t - 1)) - 1), ((2 ** t)  * m - 1), (idx_x_pre_t - 1 + (2 * m))))
            '''
            idx_x += 1
        idx_x_pre_t += (r(n, (t - 1)) - 1)
    #
    # reverse: G-round:
    idx_x = 0
    loopAmount = int(math.floor(math.log2(n)))
    for t in range(1, loopAmount, 1):
        idx_x += (r(n, t) - 1)
    idx_x_pre_t = idx_x
    for t in range(loopAmount, 0, -1):
        idx_x_pre_t -= (r(n, (t - 1)) - 1)
        for m in range(r(n, t) - 1, -1, -1):
            if((t == 1) and (((2 * m) + 1) != (n - 1))):
                targetCircuit.ccx(z[(2 * m)], op1[(2 * m) + 1], z[(2 * m) + 1])
                #print("({}, {}) | z[{}] ^= z[{}] * op1[{}]".format(t, m, (2 * m + 1), (2 * m), (2 * m + 1)))
            elif((t != 1) and ((((2 ** t) * m) + (2 ** t) - 1) != (n - 1))):
                targetCircuit.ccx(z[((2 ** t) * m) + (2 ** (t - 1)) - 1], x[(idx_x_pre_t + (2 * m))], z[((2 ** t) * m) + (2 ** t) - 1])
                #print("({}, {}) | z[{}] ^= z[{}] * x[{}]".format(t, m, (((2 ** t) * m) + (2 ** t) - 1), (((2 ** t) * m) + (2 ** (t - 1)) - 1), (idx_x_pre_t + (2 * m))))
            '''
            if((t == 1) and  (((2 * m) + 1) == (n - 1))):
                print("({}, {}) | z[{}] ^= z[{}] * op1[{}] | STRIPPED".format(t, m, (2 * m + 1), (2 * m), (2 * m + 1)))
            elif((t != 1) and ((((2 ** t) * m) + (2 ** t) - 1) == (n - 1))):
                print("({}, {}) | z[{}] ^= z[{}] * x[{}] | STRIPPED".format(t, m, (((2 ** t) * m) + (2 ** t) - 1), (((2 ** t) * m) + (2 ** (t - 1)) - 1), (idx_x_pre_t + (2 * m))))
            '''
    #
    idx_x = 0
    loopAmount = int(math.floor(math.log2(n))) - 1
    for t in range(1, loopAmount, 1):
        idx_x += (r(n, t) - 1)
    idx_x_pre_t = 0
    for t in range(loopAmount, 0, -1):
        idx_x_pre_t = idx_x - (r(n, (t - 1)) - 1)
        for m in range(1, r(n, t), 1):
            if((t == 1) and ((2 * m + 1) != (n - 1))):
                targetCircuit.ccx(op1[2 * m], op1[(2 * m) + 1], x[idx_x])
                #print("({}, {}) | x[{}] ^= op1[{}] * op1[{}]".format(t, m, idx_x, (2 * m), (2 * m + 1)))
            elif((t != 1)):
                targetCircuit.ccx(x[idx_x_pre_t - 1 + (2 * m)], x[idx_x_pre_t + (2 * m)], x[idx_x])
                #print("({}, {}) | x[{}] ^= x[{}] * x[{}]".format(t, m, idx_x, (idx_x_pre_t - 1 + (2 * m)), (idx_x_pre_t + (2 * m))))
            idx_x += 1
        idx_x = idx_x_pre_t
    #
    # 8:
    loopAmount = n - 1
    for i in range(1, loopAmount, 1):
        if is_op2_qubit:
            targetCircuit.cx(op2[i], op1[i])
        elif op2_list[i] == 1:
            targetCircuit.x(op1[i])
    #    
    # 9:
    loopAmount = n - 1
    for i in range(0, (n - 1), 1):
        if is_op2_qubit:
            targetCircuit.ccx(op2[i], op1[i], z[i])
        elif op2_list[i] == 1:
            targetCircuit.cx(op1[i], z[i])
    #
    # 10:
    loopAmount = n - 1
    for i in range(0, loopAmount, 1):
        targetCircuit.x(op1[i])
    # END
    return

def Example1():
    '''
    Example circuit 1 shows how use the Adder.
    
    '''
    # parameter:
    c_n1 = 0b00111111 # integer: n1
    c_n2 = 0b00000001 # integer: n2
    c_sizeofn = 8 #

    # elements:
    q_n1 = QuantumRegister(c_sizeofn, "n1")
    q_n2 = QuantumRegister(c_sizeofn, "n2")
    q_carry = QuantumRegister(c_sizeofn, "carry")
    q_ancilla = QuantumRegister((c_sizeofn - w(c_sizeofn) - int(math.floor(math.log2(c_sizeofn)))), "ancilla")
    c_result = ClassicalRegister(c_sizeofn, "result")
    mainCircuit = QuantumCircuit(q_n1, q_n2, q_carry, q_ancilla, c_result)

    # assignment: n1
    for i in range(0, c_sizeofn, 1):
        if(((c_n1 >> i) & 0b1) == 1):
            mainCircuit.x(q_n1[i])
    #
    # assignment: n2
    for i in range(0, c_sizeofn, 1):
        if(((c_n2 >> i) & 0b1) == 1):
            mainCircuit.x(q_n2[i])
    #
    # addition: n1 = n1 + n2
    Adder_Draper(q_n1, q_n2, q_carry, q_ancilla, mainCircuit, True, True)

    # measure:
    mainCircuit.measure(q_n1, c_result)

    # end circuit
    simulator = AerSimulator()
    job = simulator.run(mainCircuit, shots=1)
    result = job.result()


    # analyze result:
    raw = list(result.get_counts(mainCircuit).keys())[0]
    print("{} + {} = {}".format(c_n1, c_n2, int(raw, 2)))
    pass

def Example2():
    '''
    Example circuit 2 shows how use the Adder.
    
    '''
    # parameter:
    c_n1 = 0b00111111 # integer: n1
    c_n2 = 0b00000001 # integer: n2
    c_sizeofn = 8 #

    # elements:
    q_n1 = QuantumRegister(c_sizeofn, "n1")
    q_carry = QuantumRegister(c_sizeofn, "carry")
    q_ancilla = QuantumRegister((c_sizeofn - w(c_sizeofn) - int(math.floor(math.log2(c_sizeofn)))), "ancilla")
    c_result = ClassicalRegister(c_sizeofn, "result")
    mainCircuit = QuantumCircuit(q_n1, q_carry, q_ancilla, c_result)

    # assignment: n1
    for i in range(0, c_sizeofn, 1):
        if(((c_n1 >> i) & 0b1) == 1):
            mainCircuit.x(q_n1[i])
    #
    # addition: n1 = n1 + n2
    Adder_Draper(q_n1, c_n2, q_carry, q_ancilla, mainCircuit, True, False)

    # measure:
    mainCircuit.measure(q_n1, c_result)

    # end circuit
    simulator = AerSimulator()
    job = simulator.run(mainCircuit, shots=1)
    result = job.result()


    # analyze result:
    raw = list(result.get_counts(mainCircuit).keys())[0]
    print("{} + {} = {}".format(c_n1, c_n2, int(raw, 2)))
    pass

# main():
def main():
    print("Running Example circuit 1: qubit + qubit = qubit")
    Example1()
    print("Running Example circuit 2: qubit + clbit = qubit")
    Example2()
    pass

# main() entry:
if __name__ == "__main__":
    main()
    pass
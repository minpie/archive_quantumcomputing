from projectq import MainEngine
from projectq.backends import ClassicalSimulator, ResourceCounter
from projectq.ops import (
    Measure,
    H,
    All,
    Toffoli,
    CNOT,
    X,
    Swap,
)
from projectq.meta import Loop, Compute, Uncompute, Control
import math

def OperatorGet_Q(op1):
    # op1의 값을 반환하는 함수.
    '''
    * op1 비트 배치:
    LSB ... MSB
    * 반환값 비트 배치:
    MSB ... LSB
    '''
    result = 0
    for i in range(len(op1)-1, -1, -1):
        result = result << 1
        result = result + int(op1[i])
    return result 

def OperatorSet_QI(op1, op2):
    # op1의 값을 op2로 설정하는 함수.
    '''
    * op1 비트 배치:
    LSB ... MSB
    * op2 비트 배치:
    MSB ... LSB
    '''
    for i in range(0, len(op1)):
        bit = (op2 >> i) & 0b1
        if(bit == 1):
            X | op1[i]
    return

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

def Adder_Draper(op1, op2, carry, ancilla, zerorize_carry=True, is_op2_qubit=True):
    # Drapper's adder
    '''
    Adder_Draper(op1, op2, carry, ancilla, zerorize_carry, is_op2_qubit)
    Do:
    - op1 = op1 + op2
    Operands type:
    - op1: qubit
    - op2: qubit or integer
    - carry: qubit
    - ancilla: qubit
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
            Toffoli | (op1[i], op2[i], z[i]) # z[i] = op1[i] * op2[i] = carry[i+1]
        elif op2_list[i] == 1:
            CNOT | (op1[i], z[i]) # z[i] = op1[i] * op2[i] = carry[i+1]
    #
    # 2:
    loopAmount = n
    for i in range(0, loopAmount, 1):
        if is_op2_qubit:
            CNOT | (op2[i], op1[i]) # op1[i] = op1[i] ^ op2[i] = sum[i]
        elif op2_list[i] == 1:
            X | op1[i] # op1[i] = op1[i] ^ op2[i] = sum[i]
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
                Toffoli | (op1[2 * m], op1[(2 * m) + 1], x[idx_x])
                #print("({}, {}) | x[{}] ^= op1[{}] * op1[{}]".format(t, m, idx_x, (2 * m), (2 * m + 1)))
            else:
                Toffoli | (x[idx_x_pre_t - 1 + (2 * m)], x[idx_x_pre_t + (2 * m)], x[idx_x])
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
                    Toffoli | (z[(2 * m)], op1[(2 * m) + 1], z[(2 * m) + 1])
                    #print("({}, {}) | z[{}] ^= z[{}] * op1[{}]".format(t, m, (2 * m + 1), (2 * m), (2 * m + 1)))                  
                elif((t != 1) and ((((2 ** t) * m) + (2 ** t) - 1) != (n - 1))):
                    Toffoli | (z[((2 ** t) * m) + (2 ** (t - 1)) - 1], x[(idx_x_pre_t + (2 * m))], z[((2 ** t) * m) + (2 ** t) - 1])
                    #print("({}, {}) | z[{}] ^= z[{}] * x[{}]".format(t, m, (((2 ** t) * m) + (2 ** t) - 1), (((2 ** t) * m) + (2 ** (t - 1)) - 1), (idx_x_pre_t + (2 * m))))
            else:
                if(t == 1):
                    Toffoli | (z[(2 * m)], op1[(2 * m) + 1], z[(2 * m) + 1])
                    #print("({}, {}) | z[{}] ^= z[{}] * op1[{}]".format(t, m, (2 * m + 1), (2 * m), (2 * m + 1)))
                else:
                    Toffoli | (z[((2 ** t) * m) + (2 ** (t - 1)) - 1], x[(idx_x_pre_t + (2 * m))], z[((2 ** t) * m) + (2 ** t) - 1])
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
                    Toffoli | (z[2 * m - 1], op1[2 * m], z[(2 * m)])
                    #print("({}, {}) | z[{}] ^= z[{}] * op1[{}]".format(t, m, (2 * m), (2 * m - 1), (2 * m)))
                elif((t != 1) and ((((2 ** t) * m) + (2 ** (t - 1)) - 1)  != (n - 1))):
                    Toffoli | (z[(2 ** t)  * m - 1], x[(idx_x_pre_t - 1 + (2 * m))], z[((2 ** t) * m) + (2 ** (t - 1)) - 1])
                    #print("({}, {}) | z[{}] ^= z[{}] * x[{}]".format(t, m, (((2 ** t) * m) + (2 ** (t - 1)) - 1), ((2 ** t)  * m - 1), (idx_x_pre_t - 1 + (2 * m))))
            else:
                if(t == 1):
                    Toffoli | (z[2 * m - 1], op1[2 * m], z[(2 * m)])
                    #print("({}, {}) | z[{}] ^= z[{}] * op1[{}]".format(t, m, (2 * m), (2 * m - 1), (2 * m)))
                else:
                    Toffoli | (z[(2 ** t)  * m - 1], x[(idx_x_pre_t - 1 + (2 * m))], z[((2 ** t) * m) + (2 ** (t - 1)) - 1])
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
                Toffoli | (op1[2 * m], op1[(2 * m) + 1], x[idx_x])
                #print("({}, {}) | x[{}] ^= op1[{}] * op1[{}]".format(t, m, idx_x, (2 * m), (2 * m + 1)))
            else:
                Toffoli | (x[idx_x_pre_t - 1 + (2 * m)], x[idx_x_pre_t + (2 * m)], x[idx_x])
                #print("({}, {}) | x[{}] ^= x[{}] * x[{}]".format(t, m, idx_x, (idx_x_pre_t - 1 + (2 * m)), (idx_x_pre_t + (2 * m))))
            idx_x += 1
        idx_x = idx_x_pre_t
    #
    # 4:
    loopAmount = n
    for i in range(1, loopAmount, 1):
        CNOT | (z[i - 1], op1[i])
    # 5:
    loopAmount = n - 1
    for i in range(0, loopAmount, 1):
        X | (op1[i])
    #
    # 6:
    loopAmount = n - 1
    for i in range(1, loopAmount, 1):
        if is_op2_qubit:
            CNOT | (op2[i], op1[i])
        elif op2_list[i] == 1:
            X | op1[i]
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
                Toffoli | (op1[2 * m], op1[(2 * m) + 1], x[idx_x])
                #print("({}, {}) | x[{}] ^= op1[{}] * op1[{}]".format(t, m, idx_x, (2 * m), (2 * m + 1)))
            elif((t != 1)):
                Toffoli | (x[idx_x_pre_t - 1 + (2 * m)], x[idx_x_pre_t + (2 * m)], x[idx_x])
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
                Toffoli | (z[2 * m - 1], op1[2 * m], z[(2 * m)])
                #print("({}, {}) | z[{}] ^= z[{}] * op1[{}]".format(t, m, (2 * m), (2 * m - 1), (2 * m)))
            elif((t != 1) and ((((2 ** t) * m) + (2 ** (t - 1)) - 1)  != (n - 1))):
                Toffoli | (z[(2 ** t)  * m - 1], x[(idx_x_pre_t - 1 + (2 * m))], z[((2 ** t) * m) + (2 ** (t - 1)) - 1])
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
                Toffoli | (z[(2 * m)], op1[(2 * m) + 1], z[(2 * m) + 1])
                #print("({}, {}) | z[{}] ^= z[{}] * op1[{}]".format(t, m, (2 * m + 1), (2 * m), (2 * m + 1)))
            elif((t != 1) and ((((2 ** t) * m) + (2 ** t) - 1) != (n - 1))):
                Toffoli | (z[((2 ** t) * m) + (2 ** (t - 1)) - 1], x[(idx_x_pre_t + (2 * m))], z[((2 ** t) * m) + (2 ** t) - 1])
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
                Toffoli | (op1[2 * m], op1[(2 * m) + 1], x[idx_x])
                #print("({}, {}) | x[{}] ^= op1[{}] * op1[{}]".format(t, m, idx_x, (2 * m), (2 * m + 1)))
            elif((t != 1)):
                Toffoli | (x[idx_x_pre_t - 1 + (2 * m)], x[idx_x_pre_t + (2 * m)], x[idx_x])
                #print("({}, {}) | x[{}] ^= x[{}] * x[{}]".format(t, m, idx_x, (idx_x_pre_t - 1 + (2 * m)), (idx_x_pre_t + (2 * m))))
            idx_x += 1
        idx_x = idx_x_pre_t
    #
    # 8:
    loopAmount = n - 1
    for i in range(1, loopAmount, 1):
        if is_op2_qubit:
            CNOT | (op2[i], op1[i])
        elif op2_list[i] == 1:
            X | op1[i]
    #    
    # 9:
    loopAmount = n - 1
    for i in range(0, (n - 1), 1):
        if is_op2_qubit:
            Toffoli | (op2[i], op1[i], z[i])
        elif op2_list[i] == 1:
            CNOT | (op1[i], z[i])
    #
    # 10:
    loopAmount = n - 1
    for i in range(0, loopAmount, 1):
        X | (op1[i])
    # END
    return

def check(a, b):
    objEngine = MainEngine(ClassicalSimulator())


    # 수행:
    n = max(GetBitLen(a), GetBitLen(b))
    num1 = objEngine.allocate_qureg(n)
    num2 = objEngine.allocate_qureg(n)
    carry = objEngine.allocate_qureg(n)
    ancilla = objEngine.allocate_qureg(n - w(n) - int(math.floor(math.log2(n))))
    OperatorSet_QI(num1, a)
    OperatorSet_QI(num2, b)
    Adder_Draper(num1, num2, carry, ancilla)
    All(Measure) | num1
    All(Measure) | num2
    All(Measure) | carry
    All(Measure) | ancilla
    objEngine.flush()

    
    res1 = 0
    result_sum = OperatorGet_Q(num1)
    if(result_sum != ((a + b) % (2 ** n))):
        res1 += 1
    if(OperatorGet_Q(carry) != 0):
        res1 += 2
    if(OperatorGet_Q(ancilla) != 0):
        res1 += 4
    out = [res1, "({}, {}) | n = {} | sum = {} | carry = {} | ancilla = {}".format(a, b, n, result_sum, OperatorGet_Q(carry), OperatorGet_Q(ancilla))]
    return out

def main():
    cnt = 0
    n = 32
    for i in range(0, n, 1):
        for j in range(0, n, 1):
            a = (2 ** i) - 1
            b = (2 ** j) - 1
            resorg = check(a, b)
            res = resorg[0]
            out = ""
            if((res & 1) == 1):
                out += "Wrong sum!        | "
            else:
                out += "                  | "
            if((res & 2) == 2):
                out += "Non-Zero carry!   | "
            else:
                out += "                  | "
            if((res & 4) == 4):
                out += "Non-Zero ancilla! | "
            else:
                out += "                  | "
            out += resorg[1]
            if(res > 0):
                cnt += 1
                print(out)
    print("total wrong result: {} / {}".format(cnt, n * n))
    pass

#### main() entry
if __name__ == "__main__":
    main()
    pass
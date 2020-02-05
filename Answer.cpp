//------------------------------------------------------------------------------
/// @file
/// @author   ハル研究所プログラミングコンテスト実行委員会
///
/// @copyright  Copyright (c) 2019 HAL Laboratory, Inc.
/// @attention  このファイルの利用は、同梱のREADMEにある
///             利用条件に従ってください。
//------------------------------------------------------------------------------

#include "Answer.hpp"
#include "Random.hpp"
#include <vector>
#include <queue>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <iterator>
#include <random>
#include <map>
#include <fstream>
#include <string>
using std::cerr;
using std::endl;
using std::map;
using std::priority_queue;
using std::vector;

//------------------------------------------------------------------------------
namespace hpc
{
const int FINISH = 1 << 30;
Random random(RandomSeed::DefaultSeed());
//------------------------------------------------------------------------------
/// 二点のマンハッタン距離を計算する
/// @param[in] p1,p2 距離を計算する2点
/// @param[out] distance 二点間のマンハッタン距離
int manhattanDistance(Point p1, Point p2)
{
    int distance = 0;
    distance += Math::Abs(p1.x - p2.x);
    distance += Math::Abs(p1.y - p2.y);
    return distance;
}
class State
{
public:
    vector<Point> simTurtlePos;
    // vector<int> lossTurn(simTurtlePos.count(), 0);
    vector<int> lastTurn;
    vector<int> prevEatenFoodIndex;
    vector<int> nowEatenFoodIndex;
    vector<bool> wasAlreadyEaten;
    int curIndex;
    State()
    {
    }
    ~State() {}
}; // namespace hpc
const long long MOD = 1e9 + 7;
map<long long, State> stateCache;
long long vec2hash(vector<int> &v, int n)
{
    long long base = 10007;
    long long ret = 0;
    for (int i = 0; i < n; ++i)
    {
        ret += (v[i] * base) % MOD;
        ret %= MOD;
        base = (base * 10007) % MOD;
    }
    return ret;
}

class Tactics
{
public:
    vector<int> visitOrder;
    int score;
    int worstFoodIndex;
    Tactics()
    {
    }
    ~Tactics() {}
}; // namespace hpc

bool operator>(const Tactics &t1, const Tactics &t2)
{
    return t1.score > t2.score;
};

Tactics bestTactics;
vector<int> turtlesTarget;
int stageNum = -1;
int imp = 0;
int colab = 0;
vector<int> incval;
bool isSimulate = false;
int bestTurn;
int gapTurnSum = 0;
vector<vector<int>> moveOrder;
vector<vector<vector<bool>>> canVisitDuringMove;
map<std::pair<int, int>, vector<int>> foodsIndexDuringMove;
void genMoveOrder(Tactics &tact, const Stage &aStage);
vector<int> genOptimaticDuringMovement(int start, int height, int end, const Stage &aStage, vector<bool> &wasAlreadyEaten);
map<std::pair<std::pair<int, int>, int>, vector<int>> optimaticDuringMovementCache;

int evalTacticsBySimulate(Tactics &tact, const Stage &aStage, long long restoreHash, int limit);
void evalTactics(Tactics &aTact, const Stage &aStage)
{
    int ret = 0;
    int foodSize = (int)aTact.visitOrder.size();
    priority_queue<int, vector<int>, std::greater<int>> pq;
    for (Point p : aStage.turtlePositions())
    {
        pq.push(manhattanDistance(p, aStage.foods()[0].pos()));
    }
    for (int i = 1; i < aStage.foods()[0].height(); ++i)
    {
        ret += pq.top();
        pq.pop();
    }

    for (int i = 1; i < foodSize; ++i)
    {
        ret += manhattanDistance(aStage.foods()[aTact.visitOrder[i - 1]].pos(), aStage.foods()[aTact.visitOrder[i % foodSize]].pos());
    }
    //Food firstFood = aStage.foods()[aTact.visitOrder[0]];
    //ret += firstFood.height() / 5;
    //ret += abs(firstFood.pos().x - 30) / 2;
    //ret += abs(firstFood.pos().y - 15) / 2;

    aTact.score = ret;
}
//vector<int> genBaseMoveOrder(const Stage &aStage, int initFoodIndex)
//{

//}

void twoOpt(Tactics &tact, const Stage &aStage, const int loopNum)
{

    for (int i = 0; i < loopNum; ++i)
    {
        int r1, r2;
        do
        {
            r1 = random.randTerm(aStage.foods().count() - 1);
            r2 = random.randTerm(aStage.foods().count() - 1);
        } while (r1 == r2);
        if (r1 > r2)
        {
            std::swap(r1, r2);
        }
        int newScore = tact.score;
        int diff = 0;
        //diff += abs(aStage.foods()[tact.visitOrder[r1]].height() - aStage.foods()[tact.visitOrder[r1 + 1]].height());
        //diff += abs(aStage.foods()[tact.visitOrder[r2]].height() - aStage.foods()[tact.visitOrder[r2 + 1]].height());

        diff -= manhattanDistance(aStage.foods()[tact.visitOrder[r1]].pos(), aStage.foods()[tact.visitOrder[r1 + 1]].pos());
        diff -= manhattanDistance(aStage.foods()[tact.visitOrder[r2]].pos(), aStage.foods()[tact.visitOrder[r2 + 1]].pos());
        diff += manhattanDistance(aStage.foods()[tact.visitOrder[r1]].pos(), aStage.foods()[tact.visitOrder[r2]].pos());
        diff += manhattanDistance(aStage.foods()[tact.visitOrder[r2 + 1]].pos(), aStage.foods()[tact.visitOrder[r1 + 1]].pos());
        //diff -= abs(aStage.foods()[tact.visitOrder[r1]].height() - aStage.foods()[tact.visitOrder[r2 + 1]].height());
        //diff -= abs(aStage.foods()[tact.visitOrder[r2]].height() - aStage.foods()[tact.visitOrder[r1 + 1]].height());
        newScore += diff;
        //cerr << newScore << endl;
        if (newScore < tact.score)
        {

            std::reverse(tact.visitOrder.begin() + r1 + 1, tact.visitOrder.begin() + r2 + 1);

            tact.score = newScore;

            //assert(tact.score == newScore);

            // cerr << "updade! " << tact.score << endl;
        }
    }
}

Tactics searchTactics(vector<int> &targetedFoods, const Stage &aStage)
{
    Tactics retTactics;
    retTactics.score = (1 << 30);
    int retTurn = (1 << 30);
    int foodNum = targetedFoods.size();
    for (int firstFood = 0; firstFood < foodNum; firstFood++)
    {
        Tactics tact;
        tact.visitOrder.push_back(firstFood);
        vector<bool> checked(foodNum, false);
        checked[firstFood] = true;

        for (int i = 0; i < foodNum - 1; ++i)
        {
            int bestDistante = (1 << 30);
            int next = -1;
            for (int j = 0; j < foodNum; ++j)
            {
                if (checked[j] == true)
                {
                    continue;
                }
                int curDistance = manhattanDistance(aStage.foods()[tact.visitOrder.back()].pos(), aStage.foods()[j].pos());
                if (curDistance < bestDistante)
                {
                    bestDistante = curDistance;
                    next = j;
                }
            }
            tact.visitOrder.push_back(next);
            checked[next] = true;
        }
        evalTactics(tact, aStage);
        twoOpt(tact, aStage, 6000);
        tact.score = evalTacticsBySimulate(tact, aStage, -1, 1000);

        Tactics shuff;
        int bucket = (int)tact.visitOrder.size() / 2;
        for (int i = 0; i < bucket; ++i)
        {
            shuff.visitOrder.push_back(tact.visitOrder[i]);
            shuff.visitOrder.push_back(tact.visitOrder[(int)tact.visitOrder.size() - 1 - i]);
        }
        shuff.score = evalTacticsBySimulate(shuff, aStage, -1, tact.score);
        if (tact.score > shuff.score)
        {
            tact = shuff;
        }

        int tactTurn = tact.score;
        if (tactTurn < retTurn)
        //if (retTactics.score > tact.score)
        {
            retTactics = tact;
            retTurn = tactTurn;
        }
    }

    return retTactics;
}
//------------------------------------------------------------------------------
/// コンストラクタ。
/// @detail 最初のステージ開始前に実行したい処理があればここに書きます。
clock_t starttime;
Answer::Answer()
{
    starttime = clock();
}

//------------------------------------------------------------------------------
/// デストラクタ。
/// @detail 最後のステージ終了後に実行したい処理があればここに書きます。
Answer::~Answer()
{
    cerr << "total gap turn:" << gapTurnSum << endl;
}

//------------------------------------------------------------------------------
/// 各ステージ開始時に呼び出される処理。
/// @detail 各ステージに対する初期化処理が必要ならここに書きます。
/// @param aStage 現在のステージ。
int cnt;
int improve;
vector<std::pair<int, int>> improvepair;
void Answer::initialize(const Stage &aStage)
{
    improve = 0;
    improvepair.clear();
    moveOrder = vector<vector<int>>(1000, vector<int>(aStage.turtlePositions().count(), -1));
    // canVisitDutingMove[a][b][c] -> 座標aから座標bに移動するとき、座標cを訪れるか
    foodsIndexDuringMove.clear();
    optimaticDuringMovementCache.clear();
    canVisitDuringMove = vector<vector<vector<bool>>>(120, vector<vector<bool>>(120, vector<bool>(120, false)));
    const int foodNum = aStage.foods().count();
    for (int i = 0; i < foodNum; ++i)
    {
        for (int j = 0; j < foodNum; ++j)
        {
            for (int k = 0; k < foodNum; ++k)
            {
                Point posi = aStage.foods()[i].pos();
                Point posj = aStage.foods()[j].pos();
                Point posk = aStage.foods()[k].pos();
                const int minx = std::min(posi.x, posj.x);
                const int maxx = std::max(posi.x, posj.x);
                const int miny = std::min(posi.y, posj.y);
                const int maxy = std::max(posi.y, posj.y);
                canVisitDuringMove[i][j][k] = (minx <= posk.x && posk.x <= maxx) && (miny <= posk.y && posk.y <= maxy);
                if (canVisitDuringMove[i][j][k])
                {
                    foodsIndexDuringMove[{i, j}].push_back(k);
                }
            }
        }
    }

    isSimulate = true;
    int maxFoodHeight = -1;
    for (Food f : aStage.foods())
    {
        maxFoodHeight = std::max(maxFoodHeight, f.height());
    }

    ++stageNum;
    turtlesTarget = vector<int>(aStage.turtlePositions().count(), -1);
    bestTactics.score = (1 << 30);
    Stage currentStage = aStage;
    vector<int> targetedFoods;
    for (int i = 0; i < aStage.foods().count(); ++i)
    {
        targetedFoods.push_back(i);
    }
    bestTactics = searchTactics(targetedFoods, currentStage);

    //if (aStage.foods()[bestTactics.visitOrder.front()].height() > currentStage.foods()[bestTactics.visitOrder.back()].height())
    // {
    //    std::reverse(bestTactics.visitOrder.begin(), bestTactics.visitOrder.end());
    //}
    cnt = 0;

    bestTurn = evalTacticsBySimulate(bestTactics, currentStage, -1, 1000);

    //int originTurn = bestTurn;
    //evalTactics(bestTactics, currentStage);
    //int originScore = bestTactics.score;
    int lastChange = 0;

    int beamWidth = 3;
    priority_queue<Tactics, vector<Tactics>, std::greater<Tactics>> pq;
    pq.push(bestTactics);
    //double factor = (59.0 / 240.0 / 10.0) * ((stageNum % 3) - 1);
    int count = 0;
    while (double(clock() - starttime) / CLOCKS_PER_SEC < (59.9 / 240.0) * (stageNum + 1))
    {
        auto next = pq;
        while (!pq.empty())
        {
            int psocre;
            Tactics nowTact = pq.top();
            psocre = nowTact.score;
            pq.pop();
            const int loopNum = std::max(random.randTerm(6) - 1, 1);
            for (int k = 0; k < loopNum; ++k)
            {
                int r0, r1, r2;
                //            do
                //           {
                r0 = random.randTerm(13) + 1;
                r1 = 1 + random.randTerm(aStage.foods().count() - r0 - 1); //random.randTerm(aStage.foods().count());
                r2 = r1 + r0;
                //         } while ((r2 >= aStage.foods().count() || r2 < 0)); // || (aStage.foods()[r1].height() > maxFoodHeight / 2 && aStage.foods()[r2].height() > maxFoodHeight / 2));

                std::swap(nowTact.visitOrder[r1], nowTact.visitOrder[r2]);
            }
            nowTact.score = evalTacticsBySimulate(nowTact, currentStage, -1, bestTurn);
            if (psocre > nowTact.score)
            {
                improve++;
                // cerr << r0 << endl;
                //cerr << r1 << " " << r2 << endl;
                next.push(nowTact);
            }
        }
        while ((int)pq.size() <= beamWidth && !next.empty())
        {
            pq.push(next.top());
            next.pop();
        }

        lastChange++;
        cnt++;
        count++;
    }

    //evalTactics(bestTactics, aStage);
    //incval.push_back(bestTactics.score - originScore);
    bestTactics = pq.top();
    //evalTactics(bestTactics, aStage);
    //incval.push_back(bestTactics.score - originScore);
    isSimulate = false;
    genMoveOrder(bestTactics, aStage);
    //cerr << "stage " << stageNum << " count " << cnt << " change:" << changeCount << " improve:" << originTurn - bestTurn << " increase:" << bestTactics.score - originScore << endl;
} // namespace hpc

//------------------------------------------------------------------------------
/// 各ターンのカメの行動を指定する。
/// @detail 各カメの行動を指定し、aActionsの要素に追加してください。
///         aActions[i]がi番目のカメの行動になります。
///         aActionsの要素数がaStage.turtlePositions()の要素数と異なる場合、アサートに失敗します。
/// @param[in] aStage 現在のステージ。
/// @param[out] aActions 各カメの行動を指定する配列。

void Answer::setActions(const Stage &aStage, Actions &aActions)
{

    //vector<bool> isGreedy(aStage.turtlePositions().count(), false);
    //bool firstSelected = false;
    for (int i = 0; i < aStage.turtlePositions().count(); ++i)
    {

        turtlesTarget[i] = moveOrder[aStage.turn()][i];
    }
    for (int i = 0; i < aStage.turtlePositions().count(); ++i)
    {

        if (turtlesTarget[i] == -1)
        {
            aActions.add(Action_Wait);
            continue;
        }
        int target = turtlesTarget[i];
        Point turtlePosition = aStage.turtlePositions()[i];
        Point targetFoodPosition = aStage.foods()[target].pos();
        if (turtlePosition.x < targetFoodPosition.x)
        {
            aActions.add(Action_MoveRight);
        }
        else if (turtlePosition.x > targetFoodPosition.x)
        {
            aActions.add(Action_MoveLeft);
        }
        else if (turtlePosition.y < targetFoodPosition.y)
        {
            aActions.add(Action_MoveDown);
        }
        else if (turtlePosition.y > targetFoodPosition.y)
        {
            aActions.add(Action_MoveUp);
        }
        else
        {
            aActions.add(Action_Wait);
        }
    }
    //cerr << aActions.count() << endl;
} // namespace hpc

void genMoveOrder(Tactics &tact, const Stage &aStage)
{
    auto simTurtlePos = aStage.turtlePositions();
    // vector<int> lossTurn(simTurtlePos.count(), 0);
    vector<int> lastTurn(simTurtlePos.count(), 0);
    vector<int> prevEatenFoodIndex(simTurtlePos.count(), -1);
    vector<int> nowEatenFoodIndex(simTurtlePos.count(), -1);

    vector<bool> wasAlreadyEaten(aStage.foods().count(), false);
    for (int nowFoodIndex = 0; nowFoodIndex < aStage.foods().count(); nowFoodIndex++)
    {
        int i = tact.visitOrder[nowFoodIndex];
        if (wasAlreadyEaten[i])
            continue;
        priority_queue<std::pair<int, int>, vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;
        int curHeight = aStage.foods()[i].height();
        for (int j = 0; j < simTurtlePos.count(); ++j)
        {
            pq.push({(manhattanDistance(aStage.foods()[i].pos(), simTurtlePos[j]) + lastTurn[j]), j});
        }
        //prevUsed.clear();
        vector<int> moveTurtles;
        map<int, vector<int>> moveTurtlesDuring;
        int prevConsumedTurn = 0;
        map<int, int> foodIndex;   //= nowEatenFoodIndex[pq.top().second];
        map<int, int> pointHeight; //= 0;
        while (curHeight > 0)
        {
            int t = pq.top().second;
            prevConsumedTurn = std::max(prevConsumedTurn, pq.top().first);
            pq.pop();
            pointHeight[nowEatenFoodIndex[t]]++;
            moveTurtlesDuring[nowEatenFoodIndex[t]].push_back(t);

            //if (nowEatenFoodIndex[t] == closestFoodIndex)
            // {
            //    closestPointHeight++;
            //   moveTurtlesDuring.push_back(t);
            //}
            moveTurtles.push_back(t);
            //     lossTurn[t] = 0; //std::max(0, -pq.top().first);
            prevEatenFoodIndex[t] = nowEatenFoodIndex[t];
            nowEatenFoodIndex[t] = i;

            simTurtlePos[t] = aStage.foods()[i].pos();
            curHeight--;
        }
        for (auto during : pointHeight)
        {

            //cerr << during.first << " " << during.second << endl;
            int closestPointHeight = during.second;
            int closestFoodIndex = during.first;
            if (closestFoodIndex == -1)
                continue;
            auto duringMoveFoodsIndex = genOptimaticDuringMovement(closestFoodIndex, closestPointHeight, i, aStage, wasAlreadyEaten);
            //auto duringMoveFoodsIndex = optimaticDuringMovementCache[{{closestFoodIndex, i}, closestPointHeight}] = genOptimaticDuringMovement(closestFoodIndex, closestPointHeight, i, aStage, wasAlreadyEaten);
            //cerr << duringMoveFoodsIndex.size() << endl;
            int cur = closestFoodIndex;
            // cerr << (int)duringMoveFoodsIndex.size() << endl;
            for (int x : duringMoveFoodsIndex)
            {
                //assert(canVisitDuringMove[cur][i][x]);
                int dist = manhattanDistance(aStage.foods()[cur].pos(), aStage.foods()[x].pos());
                for (int t : moveTurtlesDuring[closestFoodIndex])
                {
                    for (int k = 0; k < dist; ++k)
                    {
                        moveOrder[lastTurn[t]][t] = x;
                        lastTurn[t]++;
                    }
                }
                cur = x;
                wasAlreadyEaten[x] = true;
            }
            //cerr << foodsIndexDuringMove[{closestFoodIndex, nowFoodIndex}].size() << endl;
        }
        for (int t : moveTurtles)
        {
            //assert(lastTurn[t] <= prevConsumedTurn);
            while (lastTurn[t] < prevConsumedTurn)
            {

                moveOrder[lastTurn[t]][t] = i;
                lastTurn[t]++;
            }
        }
    }
}
vector<int> genOptimaticDuringMovement(int start, int height, int end, const Stage &aStage, vector<bool> &wasAlreadyEaten)
{
    vector<int> ret;
    std::deque<vector<int>> que;
    que.push_back({start});
    int lim = foodsIndexDuringMove[{start, end}].size();
    while (!que.empty() && lim >= 0)
    {
        auto v = que.front();
        que.pop_front();
        if (ret.size() < v.size())
        {
            ret = v;
        }
        if (v.back() == end)
            continue;

        for (int x : foodsIndexDuringMove[{v.back(), end}])
        {
            if (wasAlreadyEaten[x] == false && aStage.foods()[x].height() <= height && v.back() != x)
            {
                //assert(canVisitDuringMove[v.back()][end][x]);
                auto tmp = v;
                tmp.push_back(x);
                que.push_front(tmp);
            }
        }
        lim--;
    }
    return ret;
}
// シミュレーション用
int evalTacticsBySimulate(Tactics &tact, const Stage &aStage, long long restoreHash, int limit)
{
    int ret = -1;
    auto simTurtlePos = aStage.turtlePositions();
    // vector<int> lossTurn(simTurtlePos.count(), 0);
    vector<int> lastTurn(simTurtlePos.count(), 0);
    vector<int> prevEatenFoodIndex(simTurtlePos.count(), -1);
    vector<int> nowEatenFoodIndex(simTurtlePos.count(), -1);
    vector<bool> wasAlreadyEaten(aStage.foods().count(), false);
    int startIndex = 0;

    /*if (restoreHash != -1)
    {
        if (stateCache[restoreHash].simTurtlePos.size() != 0)
        {
            lastTurn = stateCache[restoreHash].lastTurn;
            prevEatenFoodIndex = stateCache[restoreHash].prevEatenFoodIndex;
            nowEatenFoodIndex = stateCache[restoreHash].nowEatenFoodIndex;
            wasAlreadyEaten = stateCache[restoreHash].wasAlreadyEaten;
            startIndex = stateCache[restoreHash].curIndex;
        }
    }*/

    //priority_queue<std::pair<int, int>, vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> actualTime;

    for (int nowFoodIndex = startIndex; nowFoodIndex < aStage.foods().count(); nowFoodIndex++)
    {
        int i = tact.visitOrder[nowFoodIndex];
        if (wasAlreadyEaten[i])
            continue;
        priority_queue<std::pair<int, int>, vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;
        int curHeight = aStage.foods()[i].height();
        for (int j = 0; j < simTurtlePos.count(); ++j)
        {
            pq.push({(manhattanDistance(aStage.foods()[i].pos(), simTurtlePos[j]) + lastTurn[j]), j});
        }
        //prevUsed.clear();
        vector<int> moveTurtles;
        map<int, vector<int>> moveTurtlesDuring;

        map<int, int> foodIndex;   //= nowEatenFoodIndex[pq.top().second];
        map<int, int> pointHeight; //= 0;
        int prevConsumedTurn = 0;
        while (curHeight > 0)
        {
            int t = pq.top().second;
            prevConsumedTurn = std::max(prevConsumedTurn, pq.top().first);
            pq.pop();
            pointHeight[nowEatenFoodIndex[t]]++;
            moveTurtlesDuring[nowEatenFoodIndex[t]].push_back(t);

            //if (nowEatenFoodIndex[t] == closestFoodIndex)
            // {
            //    closestPointHeight++;
            //   moveTurtlesDuring.push_back(t);
            //}
            moveTurtles.push_back(t);
            //     lossTurn[t] = 0; //std::max(0, -pq.top().first);
            prevEatenFoodIndex[t] = nowEatenFoodIndex[t];
            nowEatenFoodIndex[t] = i;

            simTurtlePos[t] = aStage.foods()[i].pos();
            curHeight--;
        }

        wasAlreadyEaten[i] = true;
        //actualTime.push({prevConsumedTurn, i});
        for (auto during : pointHeight)
        {

            //cerr << during.first << " " << during.second << endl;
            int closestPointHeight = during.second;
            int closestFoodIndex = during.first;
            if (closestFoodIndex == -1)
                continue;
            auto duringMoveFoodsIndex = genOptimaticDuringMovement(closestFoodIndex, closestPointHeight, i, aStage, wasAlreadyEaten);
            //cerr << duringMoveFoodsIndex.size() << endl;
            // cerr << (int)duringMoveFoodsIndex.size() << endl;
            for (int x : duringMoveFoodsIndex)
            {
                if (wasAlreadyEaten[x] == true)
                    continue;
                wasAlreadyEaten[x] = true;
                // actualTime.push({manhattanDistance(aStage.foods()[closestFoodIndex].pos(), aStage.foods()[x].pos()) + lastTurn[moveTurtlesDuring[closestFoodIndex][0]], x});
            }
            //cerr << foodsIndexDuringMove[{closestFoodIndex, nowFoodIndex}].size() << endl;
        }

        for (int t : moveTurtles)
        {
            lastTurn[t] = prevConsumedTurn;
            ret = std::max(lastTurn[t], ret);
        }
        if (ret > limit)
        {
            ret = 1001;
            break;
            //cerr << "lim" << endl;
        }
        /*
        long long saveHash = vec2hash(tact.visitOrder, nowFoodIndex);
        stateCache[saveHash].lastTurn = lastTurn;
        stateCache[saveHash].prevEatenFoodIndex = prevEatenFoodIndex;
        stateCache[saveHash].nowEatenFoodIndex = nowEatenFoodIndex;
        stateCache[saveHash].wasAlreadyEaten = wasAlreadyEaten;
        stateCache[saveHash].curIndex = nowFoodIndex + 1;
    */
    }
    /*int worstDist = -1;
    int pTurn = 0;
    int ppTurn = 0;
    int pIndex = -1;
    while (!actualTime.empty())
    {
        auto p = actualTime.top();
        //    cerr << p.second << " " << p.first << endl;
        actualTime.pop();
 
        if (p.first - ppTurn > worstDist && pIndex != -1)
        {
            worstDist = p.first - pTurn;
            tact.worstFoodIndex = pIndex;
        }
        ppTurn = pTurn;
        pTurn = p.first;
        pIndex = p.second;
    }
    //cerr << tact.worstFoodIndex << endl;
    // << "---------------------------------" << endl;
    */
    return ret;
    //cerr << aActions.count() << endl;
} // namespace hpc

//------------------------------------------------------------------------------
/// 各ステージ終了時に呼び出される処理。
/// @detail 各ステージに対する終了処理が必要ならここに書きます。
/// @param aStage 現在のステージ。
void Answer::finalize(const Stage &aStage)
{
    cerr << "stage:" << stageNum << " cnt:" << cnt << " improve:" << improve << " guess:" << bestTurn << " actual:" << aStage.turn() << endl;
    gapTurnSum += aStage.turn() - bestTurn;

    /*
    const std::string fileName = "data/improve" + std::to_string(stageNum) + ".csv";
    std::ofstream ofs;
    ofs.open(fileName, std::ios::out);
    for (auto p : improvepair)
    {
        ofs << p.first << "," << p.second << endl;
    }
*/
    //colab -= std::min(0, aStage.turn() - bestTactics.score);
    //imp += std::max(0, aStage.turn() - bestTactics.score);
    //cerr << "diff:" << aStage.turn() - bestTactics.score << endl;
}

} // namespace hpc
  // EOF
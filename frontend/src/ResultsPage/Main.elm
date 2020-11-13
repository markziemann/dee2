module ResultsPage.Main exposing (..)

import Browser.Navigation as Nav
import Dict
import Dict.Extra as DExtra
import Maybe.Extra as MExtra
import ResultsPage.Helpers exposing (stageResultForDownload)
import ResultsPage.Types exposing (..)
import Routes
import SearchPage.Helpers exposing (delay, differentSearch)
import SearchPage.Types exposing (SearchParameters, SearchResults)
import Set
import SharedTypes exposing (PaginationOffset, RemoteData(..), WebData)
import Table


init : Nav.Key -> Maybe SearchParameters -> WebData SearchResults ->  Model
init navKey  maybeSearchParameters searchResults =
    { navKey = navKey
    , searchResults = searchResults
    , searchParameters = maybeSearchParameters
    , resultsTableQuery = ""
    , resultsTableState = Table.initialSort "id"
    , selectedResultsTableQuery = ""
    , selectedResultsTableState = Table.initialSort "id"
    , downloading = False
    , selectedResults = Dict.empty
    , resultsPendingRemoval = Set.empty
    }

gotNewSearchResults: Model -> SearchParameters ->  WebData SearchResults -> Model
gotNewSearchResults model searchParameters searchResults =
    if differentSearch model.searchParameters searchParameters then
         init model.navKey (Just searchParameters) searchResults
    else
        {model | searchParameters = Just searchParameters,  searchResults = searchResults}


onlyData : Model -> ( Model, Cmd msg)
onlyData model =
    ( model, Cmd.none)


update : Msg -> Model -> ( Model, Cmd Msg)
update msg model =
    case msg of
        ResultClicked result ->
            if Dict.member result.id model.selectedResults then
                onlyData
                    { model
                        | selectedResults = Dict.remove result.id model.selectedResults
                        , resultsPendingRemoval = Set.remove result.id model.resultsPendingRemoval
                    }

            else
                onlyData
                    { model
                        | selectedResults =
                            MExtra.unwrap model.selectedResults
                                (\value -> Dict.insert result.id value model.selectedResults)
                                (stageResultForDownload result)
                    }

        SelectedResultClicked id ->
            if Set.member id model.resultsPendingRemoval then
                onlyData { model | resultsPendingRemoval = Set.remove id model.resultsPendingRemoval }

            else
                onlyData { model | resultsPendingRemoval = Set.insert id model.resultsPendingRemoval }

        RemoveStagedSelections ->
            onlyData
                { model
                    | selectedResults =
                        DExtra.removeMany
                            (Set.intersect
                                (Set.fromList <| Dict.keys model.selectedResults)
                                model.resultsPendingRemoval
                            )
                            model.selectedResults
                    , resultsPendingRemoval = Set.empty
                }

        PageRequest searchParameters ->
            ( model, Nav.pushUrl model.navKey (Routes.searchResultsRoute searchParameters))

        SetResultsTableQuery resultsTableQuery ->
            onlyData { model | resultsTableQuery = resultsTableQuery }

        SetResultsTableState resultsTableState ->
            onlyData { model | resultsTableState = resultsTableState }

        SetSelectedResultsTableQuery selectedResultsTableQuery ->
            onlyData { model | selectedResultsTableQuery = selectedResultsTableQuery }

        SetSelectedResultsTableState selectedResultsTableState ->
            onlyData { model | selectedResultsTableState = selectedResultsTableState }

        DownloadRequested ->
            ( { model | downloading = True }, delay 5000 DownloadButtonReset)

        DownloadButtonReset ->
            onlyData { model | downloading = False }

        ShowToggleTip id ->
            --Currently this not  implemented. What needs to happen is when Toggle tip can be
            --shown (this case) we need to add flag in the table row data which will be passed
            --onto the noOverFlow column configuration which will render a button to display the
            --tool tip. The toggle tip it's self needs to be created as well. I'm thinking a
            --bootstrap dropdown. Or maybe a footer along the bottom of the page? Dunno
            onlyData model

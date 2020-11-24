module HomePage.Main exposing (..)
import Browser.Navigation as Nav
import Html exposing (Html)
import Routes
import Html exposing (..)
import Html.Attributes exposing (..)
import Html.Events exposing (onClick)

type alias Model =
    { navKey : Nav.Key }

type Msg
    = SearchRuns
    | SearchProjects

init: Nav.Key -> Model
init navKey = {navKey = navKey}

update: Msg -> Model -> (Model, Cmd msg)
update msg model =
    case msg of
        SearchRuns ->
             (model, Nav.pushUrl model.navKey (Routes.searchRunsRoute))
        SearchProjects ->
            (model, Nav.pushUrl model.navKey (Routes.searchProjectsRoute))

view: Html Msg
view =
    div [ class "card", attribute "style" "width: 18rem;" ]
    [ img [ alt "Card image cap", class "card-img-top", src "..." ]
        []
    , div [ onClick SearchRuns,class "card-body" ]
        [ h5 [ class "card-title" ]
            [ text "Card title" ]
        , p [ class "card-text" ]
            [ text "Some quick example text to build on the card title and make up the bulk of the card's content." ]
        ]
    ]